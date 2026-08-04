/*******************************************************************************
    Performance test for the warm-start solve

    Replays three flight scenarios frame by frame at 1000 Hz and solves each one
    twice: statelessly, and with a PropState attached to the rotor so that every
    element is seeded from the previous frame.

        climb      rpm and freestream both rising, -60 -> +60 m/s, up to 40k rpm
        descent    rpm and freestream both falling, +60 -> -60 m/s, from 40k rpm
        reversal   fast forward flight, then the freestream flips sign, as a
                   quick 180 would do it

    The first two sweep the whole envelope smoothly, which is the case warm
    starting is meant for. The reversal is the hard case: the operating point
    moves quickly and the freestream changes sign, so the previous frame's
    inflow angle is a poor guess and some elements change which half of the phi
    domain they solve in.

    Warm starting makes the solver path dependent (see PropState in qprop.h),
    so the two solves are not expected to agree everywhere. This test reports
    the disagreement rather than asserting it away; it only fails on a solve
    that does not converge.

    How to run:
    gcc 09_test_warm_start.c -o 09_test_warm_start -lm -Wall -Wextra
    ./09_test_warm_start [tag]           (tag defaults to "current")

    Outputs:
    09_warm_start_<tag>.csv     per-frame operating point and both solutions

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

//monotonic wall-clock with the finest resolution available, as in test 07
#if defined(_WIN32)
    #include <windows.h>
    static double now_seconds(void) {
        LARGE_INTEGER frequency, counter;
        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&counter);
        return (double)counter.QuadPart / (double)frequency.QuadPart;
    }
#elif defined(CLOCK_MONOTONIC)
    static double now_seconds(void) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (double)ts.tv_sec + 1e-9*(double)ts.tv_nsec;
    }
#else
    static double now_seconds(void) {
        return (double)clock() / (double)CLOCKS_PER_SEC;
    }
#endif

#define RATE_HZ     1000.0      //frames per second
#define MAXFRAME    2000        //longest scenario
#define NSCEN       3
#define REPEATS     5           //trajectory replays; the fastest one is reported

//smooth 0 -> 1 transition with zero slope at both ends
static double smoothstep(double x) {
    if (x <= 0.0) { return 0.0; }
    if (x >= 1.0) { return 1.0; }
    return x*x*(3.0 - 2.0*x);
}

//fill the operating point for one scenario, and return its frame count
static int build_scenario(int which, double* uinf, double* rpm, const char** name) {
    if (which == 0) {
        *name = "climb";
        int n = 2000;                                   //2.0 s
        for (int f = 0; f < n; ++f) {
            double x = f/(double)(n-1);
            uinf[f] = -60.0 + 120.0*x;                  //-60 -> +60 m/s
            rpm[f]  = 5000.0 + 35000.0*x;               //5k -> 40k rpm
        }
        return n;
    }
    if (which == 1) {
        *name = "descent";
        int n = 2000;                                   //2.0 s
        for (int f = 0; f < n; ++f) {
            double x = f/(double)(n-1);
            uinf[f] = 60.0 - 120.0*x;                   //+60 -> -60 m/s
            rpm[f]  = 40000.0 - 35000.0*x;              //40k -> 5k rpm
        }
        return n;
    }
    *name = "reversal";
    int n = 800;                                        //0.8 s
    for (int f = 0; f < n; ++f) {
        double t = f/RATE_HZ;
        //0.2 s forward, 0.4 s flipping the freestream, 0.2 s reversed
        double x = smoothstep((t - 0.2)/0.4);
        uinf[f] = 50.0 - 100.0*x;                       //+50 -> -50 m/s
        //throttle chops through the flip and comes back
        rpm[f] = 35000.0 - 10000.0*sin(PI*x);
    }
    return n;
}

static Airfoil* load_clark_y(void) {
    const char* files[12] = {
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
    return import_xfoil_polars(files, 12);
}

int main(int argc, char** argv) {
    const char* tag = (argc > 1)? argv[1] : "current";

    Airfoil* clark_y = load_clark_y();
    if (!clark_y) { printf("TEST 9 - FAILED :( could not load polars\n"); return 1; }
    Rotor* rotor = import_rotor_geometry_uiuc("../validation/hqprop_5.1x3.7x3/hqprop_5.1x3.7x3_geom.txt", clark_y, 0.12954, 3);
    if (!rotor) { printf("TEST 9 - FAILED :( could not load geometry\n"); free_airfoil(clark_y); return 1; }

    double rho = 1.225, mu = 1.81e-5, a = 340.0, tol = 1e-8;
    int itmax = 200;
    int nelems = rotor->nsections - 1;

    static double uinf[MAXFRAME], rpm[MAXFRAME];
    static double Tc[MAXFRAME], Qc[MAXFRAME], Tw[MAXFRAME], Qw[MAXFRAME];

    char csvname[128];
    snprintf(csvname, sizeof csvname, "09_warm_start_%s.csv", tag);
    FILE* csv = fopen(csvname, "w");
    if (csv) { fprintf(csv, "scenario,frame,uinf,rpm,T_cold,Q_cold,T_warm,Q_warm\n"); }

    printf("TEST 9: warm-start performance  [tag: %s]\n", tag);
    printf("  rotor hqprop 5.1x3.7x3, Clark Y, %i elements, %.0f Hz, %i replays (fastest reported)\n\n",
           nelems, RATE_HZ, REPEATS);
    printf("  %-9s %6s  %9s %9s  %8s %8s %8s   %9s  %s\n",
           "scenario", "frames", "d(Uinf)", "d(rpm)", "cold", "warm", "speed-up", "worst dT", "fails");
    printf("  %-9s %6s  %9s %9s  %8s %8s %8s   %9s\n",
           "", "", "per frame", "per frame", "ms/call", "ms/call", "", "vs cold");
    printf("  --------------------------------------------------------------------------------------------\n");

    int hard_fail = 0;

    for (int s = 0; s < NSCEN; ++s) {
        const char* name = "";
        int n = build_scenario(s, uinf, rpm, &name);

        //largest per-frame step in the operating point, which sets seed quality
        double du = 0.0, dr = 0.0;
        for (int f = 1; f < n; ++f) {
            double a1 = fabs(uinf[f]-uinf[f-1]);
            double a2 = fabs(rpm[f]-rpm[f-1]);
            if (a1 > du) { du = a1; }
            if (a2 > dr) { dr = a2; }
        }

        //---- stateless ----
        prop_state_free(rotor);
        int failc = 0;
        double best_cold = 1e300;
        for (int rep = 0; rep < REPEATS; ++rep) {
            double t0 = now_seconds();
            for (int f = 0; f < n; ++f) {
                RotorPerformance* p = qprop(rotor, uinf[f], rpm[f]*PI/30.0, tol, itmax, rho, mu, a);
                if (!p) { if (rep == 0) { ++failc; } Tc[f] = NAN; Qc[f] = NAN; continue; }
                Tc[f] = p->T;
                Qc[f] = p->Q;
                free_rotor_performance(p);
            }
            double dt = (now_seconds() - t0)/n;
            if (dt < best_cold) { best_cold = dt; }
        }

        //---- warm ----
        prop_state_new(rotor);
        int failw = 0;
        double best_warm = 1e300;
        for (int rep = 0; rep < REPEATS; ++rep) {
            prop_state_reset(rotor);        //every replay starts from cold
            double t0 = now_seconds();
            for (int f = 0; f < n; ++f) {
                RotorPerformance* p = qprop(rotor, uinf[f], rpm[f]*PI/30.0, tol, itmax, rho, mu, a);
                if (!p) { if (rep == 0) { ++failw; } Tw[f] = NAN; Qw[f] = NAN; continue; }
                Tw[f] = p->T;
                Qw[f] = p->Q;
                free_rotor_performance(p);
            }
            double dt = (now_seconds() - t0)/n;
            if (dt < best_warm) { best_warm = dt; }
        }
        prop_state_free(rotor);

        //---- agreement ----
        double worst = 0.0;
        int worst_f = 0, ndiff = 0, ncmp = 0;
        for (int f = 0; f < n; ++f) {
            if (Tc[f] != Tc[f] || Tw[f] != Tw[f]) { continue; }
            ++ncmp;
            double d = fabs(Tw[f]-Tc[f])/fmax(fabs(Tc[f]), 1e-6);
            if (d > 1e-9) { ++ndiff; }
            if (d > worst) { worst = d; worst_f = f; }
            if (csv) {
                fprintf(csv, "%s,%d,%.4f,%.1f,%.9g,%.9g,%.9g,%.9g\n",
                        name, f, uinf[f], rpm[f], Tc[f], Qc[f], Tw[f], Qw[f]);
            }
        }

        printf("  %-9s %6d  %8.3f%s %9.1f  %8.4f %8.4f %7.2fx   %8.3f%%  %d/%d\n",
               name, n, du, " m/s", dr, best_cold*1e3, best_warm*1e3,
               best_cold/best_warm, worst*100.0, failc, failw);
        printf("  %-9s %6s  worst disagreement at frame %d: Uinf %+.1f m/s, %.0f rpm, %d of %d frames differ\n",
               "", "", worst_f, uinf[worst_f], rpm[worst_f], ndiff, ncmp);

        if (failc > 0 || failw > 0) { hard_fail = 1; }
    }

    if (csv) { fclose(csv); printf("\n  wrote %s\n", csvname); }

    printf("\n  Warm starting is path dependent: near a bifurcation both halves of the\n");
    printf("  phi domain admit a root, and which one is found depends on how the\n");
    printf("  operating point was reached. The disagreement above is that effect, not\n");
    printf("  a convergence failure - every reported solve is a converged root.\n");

    free_rotor(rotor);
    free_airfoil(clark_y);

    if (hard_fail) {
        printf("\nTEST 9 - FAILED :( (a solve did not converge)\n");
        return 1;
    }
    printf("\nTEST 9 - PASSED :)\n");
    return 0;
}
