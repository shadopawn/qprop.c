/*******************************************************************************
    Testing program for the qprop() function

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

int main() {
    //load NACA-4412 polars
    const char* filenames1[10] = {
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.030_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.040_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.060_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.080_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.100_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.130_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.160_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.200_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.300_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca4412_Ncrit=6/NACA 4412_T1_Re0.500_M0.00_N6.0.txt"
    };
    Airfoil* naca4412 = import_xfoil_polars(filenames1, 10);

    //load propeller geometry from APC file
    Rotor* apc10x7sf = import_rotor_geometry_apc("../validation/apc_10x7sf/10x7SF-PERF.PE0", naca4412);
    double Omega = 6014*M_PI/30;
    double tol = 1e-6;
    int itmax = 100;
    double rho = 1.225;
    double mu = 1.81e-5;
    double a = 0.0;

    double Uinf = 1.2729633333333334;

    clock_t start, end;
    int iterations = 1000;
    
    start = clock();
    RotorPerformance* perf1;
    for (int i=0; i<iterations; ++i) {
        perf1 = qprop(apc10x7sf, Uinf, Omega, tol, itmax, rho, mu, a);
    }
    end = clock();

    double run_time_s = ((double) (end - start)) / CLOCKS_PER_SEC;
    double run_time_ms = (run_time_s / iterations) * 1000;
    
    printf("TEST 7.1:\n");
    //printf("    Uinf: %f\n", Uinf);
    //printf("    Thrust: %f\n", perf1->T);
    //printf("    Torque: %f\n", perf1->Q);

    printf("qprop ran in %f miliseconds\n", run_time_ms);

    free_rotor(apc10x7sf);
    free_airfoil(naca4412);
    free_rotor_performance(perf1);
    return 0;
}
