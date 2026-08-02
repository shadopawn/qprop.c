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
    const char* filenames[5] = {
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-50000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-100000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-200000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-500000.txt",
        "../webgui/airfoil_polars/naca2411_Ncrit=9/xf-naca2411-il-1000000.txt",
    };
    Airfoil* naca2411 = import_xfoil_polars(filenames, 5);

    //load NACA-2411 polars
    const char* filenames_360[5] = {
        "../webgui/airfoil_polars/naca2411_360_Ncrit=6/NACA_2411_NACA_2411_Re0.050_M0.00_N6.0_360_V.txt",
        "../webgui/airfoil_polars/naca2411_360_Ncrit=6/NACA_2411_NACA_2411_Re0.100_M0.00_N6.0_360_V.txt",
        "../webgui/airfoil_polars/naca2411_360_Ncrit=6/NACA_2411_NACA_2411_Re0.200_M0.00_N6.0_360_V.txt",
        "../webgui/airfoil_polars/naca2411_360_Ncrit=6/NACA_2411_NACA_2411_Re0.500_M0.00_N6.0_360_M.txt",
        "../webgui/airfoil_polars/naca2411_360_Ncrit=6/NACA_2411_NACA_2411_Re1.000_M0.00_N6.0_360_M.txt",
    };
    Airfoil* naca2411_360 = import_xfoil_polars(filenames_360, 5);

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


    const char* naca0012_filenames[12] = {
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.030_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.040_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.060_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.080_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.100_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.130_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.160_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.200_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.300_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re0.500_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re1.000_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/naca0012_Ncrit=6/NACA 0012_T1_Re3.000_M0.00_N6.0.txt"
    };
    Airfoil* naca0012 = import_xfoil_polars(naca0012_filenames, 12);

    const char* eppler_e63_filenames[12] = {
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.030_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.040_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.060_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.080_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.100_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.130_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.160_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.200_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.300_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re0.500_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re1.000_M0.00_N6.0.txt",
        "../webgui/airfoil_polars/eppler_e63_Ncrit=6/E63_T1_Re3.000_M0.00_N6.0.txt"
    };
    Airfoil* eppler_e63 = import_xfoil_polars(eppler_e63_filenames, 12);

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
    Rotor* hqprop_5137 = import_rotor_geometry_uiuc("../validation/hqprop_5.1x3.7x3/hqprop_5.1x3.7x3_geom.txt", clark_y, 5.1 / 39.37, 3);

    //load propeller geometry from APC file
    Rotor* apc5x37e3 = import_rotor_geometry_apc("../validation/apc_5x37e-3/5x37E-3-PERF.PE0", naca4412);

    double rpm = 100;
    double Omega = rpm*M_PI/30;
    double tol = 1e-6;
    int itmax = 100;
    double rho = 1.225;
    double mu = 1.81e-5;
    double a = 340.0;

    double Uinf = 0;

    // RotorPerformance* perf1;
    // perf1 = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);

    RotorPerformance* perf1;
    while (rpm <= 10000)
    {
        Omega = rpm*M_PI/30;
        perf1 = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);

        printf("%.f %f %f\n", rpm, perf1->T, perf1->Q);

        rpm += 100;
    }

    // RotorPerformance* perf1;
    // while (Uinf <= 60)
    // {
    //     perf1 = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);

    //     printf("%f %f %f\n", Uinf, perf1->T, perf1->Q);

    //     Uinf++;
    // }
    

    // clock_t start, end;
    // int iterations = 1000;
    
    // start = clock();
    // RotorPerformance* perf1;
    // for (int i=0; i<iterations; ++i) {
    //     perf1 = qprop(hqprop_5137, Uinf, Omega, tol, itmax, rho, mu, a);
    // }
    // end = clock();

    // double run_time_s = ((double) (end - start)) / CLOCKS_PER_SEC;
    // double run_time_ms = (run_time_s / iterations) * 1000;
    
    // printf("TEST 7.1:\n");
    // printf("    Uinf: %f\n", Uinf);
    // printf("    RPM: %f\n", rpm);
    // printf("    Thrust N: %f\n", perf1->T);
    // printf("    Thrust gf: %f\n", perf1->T * 101.9716212978f);
    // printf("    Torque: %f\n", perf1->Q);

    // printf("qprop ran in %f miliseconds\n", run_time_ms);

    free_rotor(hqprop_5137);
    free_airfoil(naca2411);
    free_rotor_performance(perf1);
    return 0;
}
