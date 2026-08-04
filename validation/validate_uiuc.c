/* Validate qprop against the UIUC wind-tunnel data in validation/.
   Run from the validation/ directory.

   Build against either solver via -DQPROP_SRC:
     gcc validate_uiuc.c -DQPROP_SRC='"../src/qprop.c"' -o new -lm
     gcc validate_uiuc.c -DQPROP_SRC='"<old>/qprop.c"'   -o old -lm

   Every UIUC data point is reproduced exactly: the static sets fix RPM, the
   performance sets fix J, so Uinf = J*n*D is back-computed per point. No
   interpolation of either curve is involved.

   Air properties and solver settings match the published validation scripts:
   rho=1.225, mu=1.81e-5, a=340, tol=1e-8, itmax=200, Uinf=0.01 for static. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef QPROP_SRC
#define QPROP_SRC "../src/qprop.c"
#endif
#include QPROP_SRC

#define MAXPTS 64
#define MAXFILES 8

typedef struct {
    const char* name;
    const char* geom;
    const char* polar_dir;
    const char* polar_files[16];
    int npolars;
    const char* static_file;
    const char* perf_files[MAXFILES];
    double perf_rpm[MAXFILES];
    int nperf;
} Case;

/* running error accumulator */
typedef struct {
    int n;
    double sum_dCT, sum_dCT2, max_dCT;
    double sum_dCP, sum_dCP2, max_dCP;
    double sum_absCT, sum_absCP;
    int nfail;
} Stats;

static void stats_add(Stats* s, double CTq, double CTe, double CPq, double CPe) {
    double dCT = CTq - CTe, dCP = CPq - CPe;
    s->n++;
    s->sum_dCT += dCT;  s->sum_dCT2 += dCT*dCT;
    s->sum_dCP += dCP;  s->sum_dCP2 += dCP*dCP;
    if (fabs(dCT) > s->max_dCT) { s->max_dCT = fabs(dCT); }
    if (fabs(dCP) > s->max_dCP) { s->max_dCP = fabs(dCP); }
    s->sum_absCT += fabs(CTe);
    s->sum_absCP += fabs(CPe);
}

static void stats_print(const char* label, Stats* s) {
    if (s->n == 0) { printf("  %-34s no points\n", label); return; }
    double rmsCT = sqrt(s->sum_dCT2/s->n), rmsCP = sqrt(s->sum_dCP2/s->n);
    double mCT = s->sum_absCT/s->n, mCP = s->sum_absCP/s->n;
    printf("  %-34s %3d  %8.5f %+8.5f %6.2f%%   %8.5f %+8.5f %6.2f%%",
           label, s->n,
           rmsCT, s->sum_dCT/s->n, (mCT > 1e-9)? rmsCT/mCT*100.0 : 0.0,
           rmsCP, s->sum_dCP/s->n, (mCP > 1e-9)? rmsCP/mCP*100.0 : 0.0);
    if (s->nfail) { printf("   %d FAILED", s->nfail); }
    printf("\n");
}

static void stats_merge(Stats* d, const Stats* s) {
    d->n += s->n;
    d->sum_dCT += s->sum_dCT;  d->sum_dCT2 += s->sum_dCT2;
    d->sum_dCP += s->sum_dCP;  d->sum_dCP2 += s->sum_dCP2;
    d->sum_absCT += s->sum_absCT; d->sum_absCP += s->sum_absCP;
    if (s->max_dCT > d->max_dCT) { d->max_dCT = s->max_dCT; }
    if (s->max_dCP > d->max_dCP) { d->max_dCP = s->max_dCP; }
    d->nfail += s->nfail;
}

/* read a whitespace-separated table with a one-line header */
static int read_table(const char* path, double col[][4], int ncols, int max) {
    FILE* f = fopen(path, "r");
    if (!f) { printf("  ERROR: cannot open %s\n", path); return 0; }
    char line[512];
    if (!fgets(line, sizeof line, f)) { fclose(f); return 0; }
    int n = 0;
    while (n < max && fgets(line, sizeof line, f)) {
        double v[4] = {0,0,0,0};
        int got = sscanf(line, "%lf %lf %lf %lf", &v[0], &v[1], &v[2], &v[3]);
        if (got >= ncols) {
            for (int k = 0; k < 4; ++k) { col[n][k] = v[k]; }
            ++n;
        }
    }
    fclose(f);
    return n;
}

int main(int argc, char** argv) {
    const char* tag = (argc > 1)? argv[1] : "run";

    static const char* naca4412[10] = {
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
    static const char* clarky[12] = {
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
        "../webgui/airfoil_polars/clark_y_Ncrit=7/CLARK Y AIRFOIL_T1_Re0.500_M0.00_N7.0.txt"
    };

    Case cases[3];
    memset(cases, 0, sizeof cases);

    cases[0].name = "APC 10x7SF";
    cases[0].geom = "apc_10x7sf/10x7SF-PERF.PE0";
    cases[0].npolars = 10;
    for (int i = 0; i < 10; ++i) { cases[0].polar_files[i] = naca4412[i]; }
    cases[0].static_file = "apc_10x7sf/uiuc_data/apcsf_10x7_static_kt0827.txt";
    cases[0].perf_files[0] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0828_3008.txt"; cases[0].perf_rpm[0] = 3008;
    cases[0].perf_files[1] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0829_4011.txt"; cases[0].perf_rpm[1] = 4011;
    cases[0].perf_files[2] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0830_3999.txt"; cases[0].perf_rpm[2] = 3999;
    cases[0].perf_files[3] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0831_5003.txt"; cases[0].perf_rpm[3] = 5003;
    cases[0].perf_files[4] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0832_5006.txt"; cases[0].perf_rpm[4] = 5006;
    cases[0].perf_files[5] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0833_6006.txt"; cases[0].perf_rpm[5] = 6006;
    cases[0].perf_files[6] = "apc_10x7sf/uiuc_data/apcsf_10x7_kt0834_6014.txt"; cases[0].perf_rpm[6] = 6014;
    cases[0].nperf = 7;

    cases[1].name = "APC 16x8E";
    cases[1].geom = "apc_16x8e/16x8E-PERF.PE0";
    cases[1].npolars = 10;
    for (int i = 0; i < 10; ++i) { cases[1].polar_files[i] = naca4412[i]; }
    cases[1].static_file = "apc_16x8e/uiuc_data/apce_16x8_static_2150od.txt";
    cases[1].perf_files[0] = "apc_16x8e/uiuc_data/apce_16x8_2154od_4968.txt"; cases[1].perf_rpm[0] = 4968;
    cases[1].perf_files[1] = "apc_16x8e/uiuc_data/apce_16x8_2155od_5027.txt"; cases[1].perf_rpm[1] = 5027;
    cases[1].nperf = 2;

    cases[2].name = "APC 4.2x4";
    cases[2].geom = "apc_4.2x4/42x4-PERF.PE0";
    cases[2].npolars = 12;
    for (int i = 0; i < 12; ++i) { cases[2].polar_files[i] = clarky[i]; }
    cases[2].static_file = "apc_4.2x4/uiuc_data/apcff_4.2x4_static_0615rd.txt";
    cases[2].perf_files[0] = "apc_4.2x4/uiuc_data/apcff_4.2x4_0620rd_10042.txt"; cases[2].perf_rpm[0] = 10042;
    cases[2].perf_files[1] = "apc_4.2x4/uiuc_data/apcff_4.2x4_0621rd_10071.txt"; cases[2].perf_rpm[1] = 10071;
    cases[2].nperf = 2;

    const double rho = 1.225, mu = 1.81e-5, a = 340.0, tol = 1e-8;
    const int itmax = 200;

    char csvname[128];
    snprintf(csvname, sizeof csvname, "uiuc_validation_%s.csv", tag);
    FILE* csv = fopen(csvname, "w");
    if (csv) { fprintf(csv, "case,set,rpm,J,CT_exp,CT_qprop,CP_exp,CP_qprop,eta_exp,eta_qprop\n"); }

    printf("UIUC wind-tunnel validation  [%s]\n", tag);
    printf("rho=%.3f  mu=%.2e  a=%.0f  tol=%.0e  itmax=%d\n\n", rho, mu, a, tol, itmax);
    printf("  %-34s %3s  %8s %8s %7s   %8s %8s %7s\n",
           "", "N", "rmsCT", "biasCT", "nrms", "rmsCP", "biasCP", "nrms");

    Stats grand_static, grand_perf;
    memset(&grand_static, 0, sizeof grand_static);
    memset(&grand_perf, 0, sizeof grand_perf);

    for (int c = 0; c < 3; ++c) {
        Case* cs = &cases[c];
        Airfoil* af = import_xfoil_polars(cs->polar_files, cs->npolars);
        if (!af) { printf("  %s: polar load FAILED\n", cs->name); continue; }
        Rotor* rot = import_rotor_geometry_apc(cs->geom, af);
        if (!rot) { printf("  %s: geometry load FAILED\n", cs->name); free_airfoil(af); continue; }

        printf("--- %s   D = %.4f m (%.2f in), B = %d, %d sections ---\n",
               cs->name, rot->D, rot->D/0.0254, rot->B, rot->nsections);

        Stats cst, cpf;
        memset(&cst, 0, sizeof cst);
        memset(&cpf, 0, sizeof cpf);

        /* static: UIUC gives (RPM, CT, CP) */
        if (cs->static_file) {
            static double tab[MAXPTS][4];
            int n = read_table(cs->static_file, tab, 3, MAXPTS);
            for (int i = 0; i < n; ++i) {
                double rpm = tab[i][0];
                double Omega = rpm*PI/30.0;
                RotorPerformance* p = qprop(rot, 0.01, Omega, tol, itmax, rho, mu, a);
                if (!p) { ++cst.nfail; continue; }
                stats_add(&cst, p->CT, tab[i][1], p->CP, tab[i][2]);
                if (csv) {
                    fprintf(csv, "%s,static,%.0f,0,%.6f,%.6f,%.6f,%.6f,,\n",
                            cs->name, rpm, tab[i][1], p->CT, tab[i][2], p->CP);
                }
                free_rotor_performance(p);
            }
            stats_print("static (RPM sweep)", &cst);
        }

        /* performance: UIUC gives (J, CT, CP, eta) at a fixed RPM per file */
        for (int k = 0; k < cs->nperf; ++k) {
            static double tab[MAXPTS][4];
            int n = read_table(cs->perf_files[k], tab, 4, MAXPTS);
            double rpm = cs->perf_rpm[k];
            double nrev = rpm/60.0;
            double Omega = rpm*PI/30.0;
            Stats fs;
            memset(&fs, 0, sizeof fs);
            for (int i = 0; i < n; ++i) {
                double J = tab[i][0];
                double Uinf = J*nrev*rot->D;
                RotorPerformance* p = qprop(rot, Uinf, Omega, tol, itmax, rho, mu, a);
                if (!p) { ++fs.nfail; continue; }
                stats_add(&fs, p->CT, tab[i][1], p->CP, tab[i][2]);
                if (csv) {
                    double eq = (fabs(p->CP) > 1e-9)? J*p->CT/p->CP : 0.0;
                    fprintf(csv, "%s,perf,%.0f,%.4f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n",
                            cs->name, rpm, J, tab[i][1], p->CT, tab[i][2], p->CP, tab[i][3], eq);
                }
                free_rotor_performance(p);
            }
            char lbl[64];
            snprintf(lbl, sizeof lbl, "perf %.0f rpm", rpm);
            stats_print(lbl, &fs);
            stats_merge(&cpf, &fs);
        }
        if (cs->nperf > 1) { stats_print("perf (all RPM)", &cpf); }
        printf("\n");

        stats_merge(&grand_static, &cst);
        stats_merge(&grand_perf, &cpf);

        free_rotor(rot);
        free_airfoil(af);
    }

    printf("--- ALL CASES ---\n");
    stats_print("static", &grand_static);
    stats_print("performance", &grand_perf);
    Stats all;
    memset(&all, 0, sizeof all);
    stats_merge(&all, &grand_static);
    stats_merge(&all, &grand_perf);
    stats_print("combined", &all);

    if (csv) { fclose(csv); printf("\nwrote %s\n", csvname); }
    return (all.nfail > 0)? 1 : 0;
}
