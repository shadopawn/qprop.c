/*******************************************************************************
    Testing program for the qprop() function

    Benchmarks qprop() over a grid of operating points, sweeping the freestream
    velocity from -60 m/s (descent / reverse inflow) to +60 m/s and the rotor
    speed up to 40000 rpm.

    Timing notes:
    - each qprop() call is timed individually with the best monotonic clock
      available on the platform, so the reported cost excludes the airfoil and
      geometry import (done once) and the free_rotor_performance() teardown
    - every operating point is warmed up before being timed, and the minimum
      over many samples is reported as the estimate of the intrinsic cost:
      the minimum is the sample least contaminated by scheduler preemption.
      The mean and the maximum are printed alongside it to expose the spread
    - the results of every call are accumulated into a volatile sink so that
      the optimizer cannot eliminate the work being measured
    - the speed of sound is set to a finite value, so the high-rpm samples
      exercise the compressibility correction in interpolate_airfoil_polars()

    Every operating point is also written to a CSV so that runs can be compared
    across a change to the solver. If a baseline is present, the run prints a
    comparison of both the results (T, Q) and the cost (time per call).

    How to run:
    gcc 07_test_performance.c -o 07_test_performance -lm -Wall -Wextra
    ./07_test_performance [tag]          (tag defaults to "current")

    Outputs:
    07_performance_<tag>.csv    per-point results and timings

    Suggested workflow for a solver change:
    ./07_test_performance before         <- record the current behaviour
    ...change the solver...
    ./07_test_performance after          <- prints deltas against the baseline

    Timings are only comparable between runs on the same machine, built with
    the same flags, and ideally with the machine otherwise idle.

    Author: Andrea Pavan
    License: MIT
*******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../src/qprop.c"

//monotonic wall-clock with the finest resolution available
//clock() is only a fallback: its tick is 1 ms on several platforms, which is
//far too coarse to time a single qprop() call
#if defined(_WIN32)
    #include <windows.h>
    static double now_seconds(void) {
        LARGE_INTEGER frequency, counter;
        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&counter);
        return (double)counter.QuadPart / (double)frequency.QuadPart;
    }
    static const char* TIMER_NAME = "QueryPerformanceCounter";
#elif defined(CLOCK_MONOTONIC)
    static double now_seconds(void) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (double)ts.tv_sec + 1e-9*(double)ts.tv_nsec;
    }
    static const char* TIMER_NAME = "clock_gettime(CLOCK_MONOTONIC)";
#else
    static double now_seconds(void) {
        return (double)clock() / (double)CLOCKS_PER_SEC;
    }
    static const char* TIMER_NAME = "clock()";
#endif

//measure the resolution of the timer, by spinning until it reports a new value
static double timer_tick_seconds(void) {
    double t0 = now_seconds();
    double t1 = now_seconds();
    while (t1 <= t0) {
        t1 = now_seconds();
    }
    return t1 - t0;
}

//measure the cost of reading the timer, so that it can be compared against the
//duration of what is being timed
static double timer_overhead_seconds(void) {
    const int n = 10000;
    double best = 1e300;
    for (int rep=0; rep<5; ++rep) {
        double t0 = now_seconds();
        for (int i=0; i<n; ++i) {
            now_seconds();
        }
        double dt = (now_seconds() - t0) / n;
        if (dt < best) {
            best = dt;
        }
    }
    return best;
}

//benchmark grid
#define UINF_MIN    (-60.0)     //minimum freestream velocity (m/s)
#define UINF_MAX    (60.0)      //maximum freestream velocity (m/s)
#define NUINF       13          //number of freestream velocity samples
#define NRPM        5           //number of rotor speed samples
#define WARMUP      10          //untimed calls per operating point
#define SAMPLES     300         //timed calls per operating point
#define NPOINTS     (NUINF*NRPM)

//below these, two runs are treated as having produced the same answer;
//both sit comfortably above the CSV's 12-significant-digit round-trip error
#define TOL_T_SAME  1e-9        //N
#define TOL_Q_SAME  1e-12       //N-m

//one grid point, kept so the sweep can be written out and diffed later
typedef struct {
    double rpm, uinf, J, T, Q;
    double tmin, tmean, tmax;   //seconds per call
    int converged;
} Row;

static void write_csv(const char* path, Row* rows, int n) {
    FILE* f = fopen(path, "w");
    if (!f) {
        printf("    WARNING: could not write %s\n", path);
        return;
    }
    //T, Q and J carry 12 significant digits so that a round-trip through the
    //file stays well below the thresholds compare() uses to call two runs
    //identical; %.8f would lose ~1e-8 on a thrust of order 10 N and make an
    //unchanged solver look like it moved every point
    fprintf(f, "rpm,uinf,J,T,Q,t_min_ms,t_mean_ms,t_max_ms,converged\n");
    for (int i=0; i<n; ++i) {
        fprintf(f, "%.0f,%.4f,%.12g,%.12g,%.12g,%.6f,%.6f,%.6f,%d\n",
                rows[i].rpm, rows[i].uinf, rows[i].J, rows[i].T, rows[i].Q,
                rows[i].tmin*1e3, rows[i].tmean*1e3, rows[i].tmax*1e3, rows[i].converged);
    }
    fclose(f);
}

static int read_csv(const char* path, Row* rows, int max) {
    FILE* f = fopen(path, "r");
    if (!f) {
        return 0;
    }
    char line[512];
    if (!fgets(line, sizeof line, f)) {      //header
        fclose(f);
        return 0;
    }
    int n = 0;
    while (n < max && fgets(line, sizeof line, f)) {
        Row* r = &rows[n];
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
                   &r->rpm, &r->uinf, &r->J, &r->T, &r->Q,
                   &r->tmin, &r->tmean, &r->tmax, &r->converged) == 9) {
            //the CSV stores milliseconds; keep seconds internally
            r->tmin  *= 1e-3;
            r->tmean *= 1e-3;
            r->tmax  *= 1e-3;
            ++n;
        }
    }
    fclose(f);
    return n;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y)? -1 : ((x > y)? 1 : 0);
}

//compare a run against a baseline: results first, then cost
static void compare(Row* cur, int ncur, Row* base, int nbase, const char* path) {
    printf("\n--- COMPARISON vs %s ---\n", path);

    double worst_dT = 0.0, worst_dT_rel = 0.0, worst_dT_u = 0.0, worst_dT_rpm = 0.0;
    double worst_dQ = 0.0, worst_dQ_rel = 0.0, worst_dQ_u = 0.0, worst_dQ_rpm = 0.0;
    int matched = 0, differing = 0, conv_changed = 0;
    static double ratios[NPOINTS];
    int nratio = 0;

    for (int i=0; i<ncur; ++i) {
        Row* b = NULL;
        for (int j=0; j<nbase; ++j) {
            if (fabs(base[j].rpm - cur[i].rpm) < 1e-6 && fabs(base[j].uinf - cur[i].uinf) < 1e-6) {
                b = &base[j];
                break;
            }
        }
        if (!b) {
            continue;
        }
        ++matched;
        if (b->converged != cur[i].converged) {
            ++conv_changed;
            printf("    convergence changed at Uinf = %+.1f m/s, %.0f rpm: %s -> %s\n",
                   cur[i].uinf, cur[i].rpm,
                   b->converged? "converged" : "failed",
                   cur[i].converged? "converged" : "failed");
        }
        if (!b->converged || !cur[i].converged) {
            continue;
        }

        double dT = cur[i].T - b->T;
        double dQ = cur[i].Q - b->Q;
        double rT = (fabs(b->T) > 1e-9)? fabs(dT)/fabs(b->T) : 0.0;
        double rQ = (fabs(b->Q) > 1e-9)? fabs(dQ)/fabs(b->Q) : 0.0;
        if (fabs(dT) > TOL_T_SAME || fabs(dQ) > TOL_Q_SAME) {
            ++differing;
        }
        if (fabs(dT) > fabs(worst_dT)) {
            worst_dT = dT; worst_dT_rel = rT; worst_dT_u = cur[i].uinf; worst_dT_rpm = cur[i].rpm;
        }
        if (fabs(dQ) > fabs(worst_dQ)) {
            worst_dQ = dQ; worst_dQ_rel = rQ; worst_dQ_u = cur[i].uinf; worst_dQ_rpm = cur[i].rpm;
        }
        if (b->tmin > 0.0 && nratio < NPOINTS) {
            ratios[nratio++] = cur[i].tmin / b->tmin;
        }
    }

    printf("    matched %i of %i grid points\n", matched, ncur);
    if (conv_changed == 0) {
        printf("    convergence: unchanged at every matched point\n");
    }

    printf("\n    results:\n");
    if (differing == 0) {
        printf("      identical at all %i points (|dT| <= %.0e N, |dQ| <= %.0e N-m)\n",
               matched, TOL_T_SAME, TOL_Q_SAME);
    } else {
        printf("      %i of %i points differ\n", differing, matched);
        printf("      largest dT = %+.6e N   (%.4f%%)  at Uinf = %+.1f m/s, %.0f rpm\n",
               worst_dT, worst_dT_rel*100.0, worst_dT_u, worst_dT_rpm);
        printf("      largest dQ = %+.6e N-m (%.4f%%)  at Uinf = %+.1f m/s, %.0f rpm\n",
               worst_dQ, worst_dQ_rel*100.0, worst_dQ_u, worst_dQ_rpm);
    }

    printf("\n    cost (min time per call, this run / baseline):\n");
    if (nratio == 0) {
        printf("      no comparable timings\n");
        return;
    }
    qsort(ratios, nratio, sizeof(double), cmp_double);
    double median = ratios[nratio/2];
    int faster = 0, slower = 0;
    for (int i=0; i<nratio; ++i) {
        if (ratios[i] < 0.98) { ++faster; }
        else if (ratios[i] > 1.02) { ++slower; }
    }
    printf("      median  %.3fx   best %.3fx   worst %.3fx\n", median, ratios[0], ratios[nratio-1]);
    printf("      faster at %i points, slower at %i, within +-2%% at %i\n",
           faster, slower, nratio - faster - slower);
    printf("      (timings only comparable on the same machine and build flags)\n");
}

int main(int argc, char** argv) {
    const char* tag = (argc > 1)? argv[1] : "current";
    const char* clark_y_filenames[12] = {
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.010_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.020_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.030_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.040_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.060_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.080_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.100_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.130_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.160_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.200_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.300_M0.00_N7.0.txt",
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.500_M0.00_N7.0.txt",
    };
    Airfoil* clark_y = import_xfoil_polars(clark_y_filenames, 12);

    //load propeller geometry from uiuc geometry file
    double D = 0.12954;
    int B = 3;
    Rotor* hqprop_5137 = import_rotor_geometry_uiuc("../validation/hqprop_5.1x3.7x3/hqprop_5.1x3.7x3_geom.txt", clark_y, D, B);
    double tol = 1e-6;
    int itmax = 100;
    double rho = 1.225;
    double mu = 1.81e-5;
    double a = 340.0;

    double rpm_samples[NRPM] = {5000.0, 10000.0, 20000.0, 30000.0, 40000.0};
    double uinf_step = (UINF_MAX - UINF_MIN) / (NUINF - 1);

    //keep the optimizer from discarding the calls being timed
    volatile double sink = 0.0;

    double tick = timer_tick_seconds();
    double overhead = timer_overhead_seconds();

    printf("TEST 7.1: qprop() benchmark\n");
    printf("    timer:     %s (%.1f ns tick, %.1f ns per reading)\n", TIMER_NAME, tick*1e9, overhead*1e9);
    printf("    rotor:     hqprop 5.1x3.7x3, D = %.5f m, B = %i, %i elements\n", D, B, hqprop_5137->nsections);
    printf("    airfoil:   Clark Y, %i polars, Ncrit = 7\n", clark_y->size);
    printf("    air:       rho = %.3f kg/m3, mu = %.3e Pa-s, a = %.1f m/s\n", rho, mu, a);
    printf("    solver:    tol = %.0e, itmax = %i\n", tol, itmax);
    printf("    sampling:  %i warm-up + %i timed calls per point, %i points\n", WARMUP, SAMPLES, NUINF*NRPM);

    //aggregate statistics over the whole grid
    double grid_min = 1e300, grid_max = 0.0, grid_mean_sum = 0.0;
    double grid_min_uinf = 0.0, grid_min_rpm = 0.0;
    double grid_max_uinf = 0.0, grid_max_rpm = 0.0;
    int npoints = 0;
    int nfailed = 0;
    static Row rows[NPOINTS];
    int nrows = 0;
    double bench_start = now_seconds();

    for (int k=0; k<NRPM; ++k) {
        double rpm = rpm_samples[k];
        double Omega = rpm*M_PI/30;
        double Vtip = Omega*D/2;

        printf("\n--- %.0f rpm  (Omega = %.1f rad/s, tip speed = %.1f m/s, tip Mach = %.2f) ---\n", rpm, Omega, Vtip, Vtip/a);
        printf("      Uinf        J       T (N)     Q (N-m)       min      mean       max\n");
        printf("     (m/s)                                       (ms)      (ms)      (ms)\n");
        printf("    ------  -------  ----------  ----------  --------  --------  --------\n");

        for (int i=0; i<NUINF; ++i) {
            double Uinf = UINF_MIN + i*uinf_step;

            //warm up the caches and check that this operating point converges
            double T = 0.0, Q = 0.0, J = 0.0;
            int converged = 1;
            for (int w=0; w<WARMUP; ++w) {
                RotorPerformance* perf = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);
                if (!perf) {
                    converged = 0;
                    break;
                }
                T = perf->T;
                Q = perf->Q;
                J = perf->J;
                free_rotor_performance(perf);
            }
            if (!converged) {
                printf("    %6.1f  %7s  %10s  %10s  %8s  %8s  %8s\n", Uinf, "-", "no conv.", "-", "-", "-", "-");
                Row* row = &rows[nrows++];
                memset(row, 0, sizeof *row);
                row->rpm = rpm;
                row->uinf = Uinf;
                row->converged = 0;
                ++nfailed;
                continue;
            }

            //time each call individually, freeing the result outside the timed region
            double tmin = 1e300, tmax = 0.0, tsum = 0.0;
            for (int s=0; s<SAMPLES; ++s) {
                double t0 = now_seconds();
                RotorPerformance* perf = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);
                double t1 = now_seconds();
                if (!perf) {
                    converged = 0;
                    break;
                }
                sink += perf->T;
                free_rotor_performance(perf);

                double dt = t1 - t0;
                if (dt < tmin) {
                    tmin = dt;
                }
                if (dt > tmax) {
                    tmax = dt;
                }
                tsum += dt;
            }
            if (!converged) {
                printf("    %6.1f  %7s  %10s  %10s  %8s  %8s  %8s\n", Uinf, "-", "no conv.", "-", "-", "-", "-");
                Row* row = &rows[nrows++];
                memset(row, 0, sizeof *row);
                row->rpm = rpm;
                row->uinf = Uinf;
                row->converged = 0;
                ++nfailed;
                continue;
            }
            double tmean = tsum / SAMPLES;

            printf("    %6.1f  %7.3f  %10.4f  %10.6f  %8.4f  %8.4f  %8.4f\n",
                   Uinf, J, T, Q, tmin*1e3, tmean*1e3, tmax*1e3);

            Row* row = &rows[nrows++];
            row->rpm = rpm;
            row->uinf = Uinf;
            row->J = J;
            row->T = T;
            row->Q = Q;
            row->tmin = tmin;
            row->tmean = tmean;
            row->tmax = tmax;
            row->converged = 1;

            if (tmin < grid_min) {
                grid_min = tmin;
                grid_min_uinf = Uinf;
                grid_min_rpm = rpm;
            }
            if (tmin > grid_max) {
                grid_max = tmin;
                grid_max_uinf = Uinf;
                grid_max_rpm = rpm;
            }
            grid_mean_sum += tmean;
            ++npoints;
        }
    }

    double bench_time = now_seconds() - bench_start;

    printf("\n--- SUMMARY ---\n");
    if (npoints > 0) {
        printf("    fastest point:   %8.4f ms/call  (Uinf = %.1f m/s, %.0f rpm)\n", grid_min*1e3, grid_min_uinf, grid_min_rpm);
        printf("    slowest point:   %8.4f ms/call  (Uinf = %.1f m/s, %.0f rpm)\n", grid_max*1e3, grid_max_uinf, grid_max_rpm);
        printf("    grid average:    %8.4f ms/call\n", (grid_mean_sum/npoints)*1e3);
    }
    printf("    operating points: %i converged, %i failed\n", npoints, nfailed);
    printf("    qprop calls:      %i\n", npoints*(WARMUP+SAMPLES));
    printf("    benchmark ran in %.2f s\n", bench_time);
    //referenced so that the accumulator is not optimized away
    printf("    (checksum: %.6e)\n", sink);

    //---- write this run out, and diff it against a baseline if one exists ----
    char csv[256];
    snprintf(csv, sizeof csv, "07_performance_%s.csv", tag);
    write_csv(csv, rows, nrows);
    printf("\n    wrote %s\n", csv);

    if (strcmp(tag, "before") == 0) {
        printf("    (this run is the baseline)\n");
    } else {
        static Row base[NPOINTS];
        int nbase = read_csv("07_performance_before.csv", base, NPOINTS);
        if (nbase > 0) {
            compare(rows, nrows, base, nbase, "07_performance_before.csv");
        } else {
            printf("    (no 07_performance_before.csv found; run with tag 'before' first)\n");
        }
    }

    free_rotor(hqprop_5137);
    free_airfoil(clark_y);

    if (nfailed > 0) {
        printf("\nTEST 7.1 - FAILED :(\n");
        return 1;
    }
    printf("\nTEST 7.1 - PASSED :)\n");
    return 0;
}
