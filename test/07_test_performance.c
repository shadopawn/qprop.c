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

    How to run:
    gcc 07_test_performance.c -o 07_test_performance -lm -Wall -Wextra
    ./07_test_performance

    Author: Andrea Pavan
    License: MIT
*******************************************************************************/
#include <math.h>
#include <stdio.h>
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

int main() {
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
                ++nfailed;
                continue;
            }
            double tmean = tsum / SAMPLES;

            printf("    %6.1f  %7.3f  %10.4f  %10.6f  %8.4f  %8.4f  %8.4f\n",
                   Uinf, J, T, Q, tmin*1e3, tmean*1e3, tmax*1e3);

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

    free_rotor(hqprop_5137);
    free_airfoil(clark_y);

    if (nfailed > 0) {
        printf("\nTEST 7.1 - FAILED :(\n");
        return 1;
    }
    printf("\nTEST 7.1 - PASSED :)\n");
    return 0;
}
