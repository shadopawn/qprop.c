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
    //load NACA-2411 polars
    const char* filenames1[5] = {
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-50000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-100000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-200000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-500000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-1000000.txt",
    };
    Airfoil* naca2411 = import_xfoil_polars(filenames1, 5);

    //load propeller geometry from uiuc geometry file
    Rotor* hqprop_5137 = import_rotor_geometry_uiuc("../validation/hqprop_5.1x3.7x3/hqprop_5.1x3.7x3_geom.txt", naca2411, 0.12954, 3);
    double rpm = 35000;
    double Omega = rpm*M_PI/30;
    double tol = 1e-6;
    int itmax = 100;
    double rho = 1.225;
    double mu = 1.81e-5;
    double a = 340.0;

    double Uinf = 42;

    clock_t start, end;
    int iterations = 1000;
    
    start = clock();
    RotorPerformance* perf1;
    for (int i=0; i<iterations; ++i) {
        perf1 = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);
    }
    end = clock();

    double run_time_s = ((double) (end - start)) / CLOCKS_PER_SEC;
    double run_time_ms = (run_time_s / iterations) * 1000;
    
    printf("TEST 7.1:\n");
    printf("    Uinf: %f\n", Uinf);
    printf("    RPM: %f\n", rpm);
    printf("    Thrust: %f\n", perf1->T);
    printf("    Torque: %f\n", perf1->Q);

    printf("qprop ran in %f miliseconds\n", run_time_ms);

    free_rotor(hqprop_5137);
    free_airfoil(naca2411);
    free_rotor_performance(perf1);
    return 0;
}
