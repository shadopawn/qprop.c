/*******************************************************************************
    Flow-state verification harness for qprop()

    Sweeps the rotor from fast descent through hover to high advance ratio and
    checks the solution against physics that any correct solver must satisfy,
    in all three regimes:

        descent / windmill brake   Uinf < 0
        vortex ring state          -2 < V/v_h < 0
        propeller branch           0 < J < J(T=0)
        turbine / windmill brake   J > J(T=0)

    This is a before/after harness, not a pass-or-die unit test. Run it against
    the current solver to record a baseline, change the solver, run it again,
    and compare. It writes a CSV and a self-contained HTML plot each run; if a
    baseline CSV is present the plot overlays both.

    How to run:
    gcc 08_test_flow_states.c -o 08_test_flow_states -lm -Wall -Wextra
    ./08_test_flow_states [tag]          (tag defaults to "current")

    Outputs:
    08_flow_states_<tag>.csv    machine-readable sweep + diagnostics
    08_flow_states_<tag>.html   plots, overlaid on 08_flow_states_before.csv

    Suggested workflow for the Ning solver swap:
    ./08_test_flow_states before         <- record the current behaviour
    ...change residual()/solver...
    ./08_test_flow_states after          <- plot overlays before vs after
    STRICT=1 ./08_test_flow_states after <- now make the physics checks binding

    Exit status:
      hard failures (non-convergence, NaN) always exit 1
      physics-admissibility failures exit 1 only when STRICT=1 is set, so the
      harness can live in run_tests.sh while the solver is still being worked on

    Author: Andrea Pavan
    License: MIT
*******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/qprop.c"

#define SWEEP_MIN   (-60.0)     //most negative freestream (m/s)
#define SWEEP_MAX   ( 60.0)     //most positive freestream (m/s)
#define SWEEP_N     241         //sweep points (0.5 m/s steps)
#define SWEEP_RPM   20000.0     //rotor speed for the swept curve

//tolerances
#define TOL_SMOOTH   0.35       //max |dT| between adjacent points, N
#define U_DEADBAND   1.0        //|Uinf| below this: sign of Uinf is meaningless
#define MOM_FLOOR    0.05       //momentum error is scaled by max(|T|, this*T_static)

//admissibility verdicts
#define ADM_NA       2          //check does not apply at this operating point
#define ADM_OK       1
#define ADM_VIOLATED 0

//a point on the sweep plus everything the checks need
typedef struct {
    double uinf, J, T, Q, P, CT, CP;
    double va;          //area-weighted induced axial velocity (m/s)
    double Wa;          //area-weighted axial velocity at the disk (m/s)
    double a;           //Ning axial induction factor, 1 - Wa/Uinf
    double vh;          //sqrt(|T|/(2 rho A))
    double Vvh;         //Uinf/vh
    double mom_abs;     //|T_be - T_annulus| in N, normalised later
    double mom_err;     //mom_abs / max(|T|, MOM_FLOOR*T_static)
    int    converged;
    int    admissible;  //ADM_* : does a steady streamtube exist (WBS region only)
    int    ning_region; //0 momentum (a<0.4), 1 empirical (0.4<a<1), 2 propeller brake (a>1)
} Pt;

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

//qprop() gives disk totals but not the induced-velocity distribution, so the
//sweep re-walks the elements to collect the diagnostics the checks need.
//It uses the library's own residual() and brent(), so it measures the solver
//under test rather than a reimplementation of it.
static void solve_point(Pt* p, Rotor* rot, double Uinf, double Omega,
                        double rho, double mu, double a_sound, double tol, int itmax) {
    int nelems = rot->nsections - 1;
    int B = rot->B;
    double R = rot->D/2;
    double A = PI*R*R;
    double n = Omega/(2*PI);

    memset(p, 0, sizeof *p);
    p->uinf = Uinf;
    p->converged = 1;

    double T = 0.0, Q = 0.0, Tann = 0.0, vaA = 0.0, Asum = 0.0;
    for (int i = 0; i < nelems; ++i) {
        Element el;
        el.c    = 0.5*(rot->sections[i].c + rot->sections[i+1].c);
        el.beta = 0.5*(rot->sections[i].beta + rot->sections[i+1].beta);
        el.r    = 0.5*(rot->sections[i].r + rot->sections[i+1].r);
        el.dr   = rot->sections[i+1].r - rot->sections[i].r;
        el.airfoil = &(rot->sections[i+1].airfoil);

        //mirror the solve inside qprop(): brake branch first for reversed
        //freestream, momentum/empirical branch otherwise or as fallback
        ResidualArgs args = {Uinf, Omega*el.r, R, B, &el, rho, mu, a_sound, 0};
        const double phi_eps = 1e-6;
        double phi = 0.0;
        int ok = 0;
        if (Uinf < 0.0) {
            args.brake = 1;
            ok = solve_bracket(&phi, -PI/2 + phi_eps, -phi_eps, tol, itmax, &args);
        }
        if (!ok) {
            args.brake = 0;
            ok = solve_bracket(&phi, phi_eps, PI/2, tol, itmax, &args);
        }
        if (!ok) { p->converged = 0; return; }
        ResidualOutput o;
        residual(&o, phi, &args);
        if (fabs(o.residual) > tol || o.W != o.W) { p->converged = 0; return; }

        double Wa_el = o.W*sin(o.phi);
        //Prandtl tip loss, same |sin phi| form the solver uses
        double s_abs = fabs(sin(o.phi));
        double F = 1.0;
        if (s_abs > 1e-9) {
            double ftip = 0.5*B*(R - el.r)/(el.r*s_abs);
            if (ftip < 0.0) { ftip = 0.0; }
            F = acos(exp(-ftip))*2.0/PI;
        }

        T += 0.5*rho*o.W*o.W*o.Cn*el.c*B*el.dr;
        Q += 0.5*rho*o.W*o.W*o.Ct*el.c*el.r*B*el.dr;

        //branch-aware annulus momentum: dT/dr = 4*pi*rho*r*|Wa|*va*F
        //|Wa| unifies CT=4a(1-a)F (momentum) and CT=4a(a-1)F (propeller brake)
        Tann += 4.0*PI*rho*el.r*fabs(Wa_el)*o.va*F*el.dr;

        double dA = 2.0*PI*el.r*el.dr;
        vaA  += o.va*dA;
        Asum += dA;
    }

    p->T = T;
    p->Q = Q;
    p->P = Q*Omega;
    p->CT = T/(rho*n*n*pow(rot->D,4));
    p->CP = 2*PI*Q/(rho*n*n*pow(rot->D,5));
    p->J  = Uinf/(n*rot->D);
    p->va = vaA/Asum;
    p->Wa = Uinf + p->va;
    p->mom_abs = fabs(T - Tann);
    p->vh = sqrt(fabs(T)/(2.0*rho*A));
    p->Vvh = (p->vh > 1e-9)? Uinf/p->vh : 0.0;

    //Ning's axial induction factor and region. Undefined at Uinf = 0.
    if (fabs(Uinf) > 1e-6) {
        p->a = 1.0 - p->Wa/Uinf;
        p->ning_region = (p->a < 0.4)? 0 : ((p->a < 1.0)? 1 : 2);
    } else {
        p->a = NAN;
        p->ning_region = 0;
    }

    //Streamtube admissibility, asserted ONLY below V/v_h = -2.
    //
    //In slow descent a rotor legitimately has Wa > 0 with Uinf < 0: the induced
    //velocity dominates and the disk still draws air through in the thrust
    //direction. That is the normal working state, not a violation.
    //
    //Below V/v_h = -2 the windmill brake state applies, the flow passes through
    //the disk in the freestream direction, and any correct solution must have
    //sign(Wa) == sign(Uinf). Between -2 and 0 lies the vortex ring state, where
    //momentum theory has no valid solution at all and nothing can be asserted.
    //The -2 threshold is exact only for a uniform actuator disk; a real blade
    //crosses branch-by-branch over a narrow band, so the assertion starts at
    //-2.1 to keep the boundary elements out of it.
    if (Uinf < -U_DEADBAND && p->Vvh <= -2.1) {
        p->admissible = (p->Wa*Uinf > 0)? ADM_OK : ADM_VIOLATED;
    } else {
        p->admissible = ADM_NA;
    }
}

static const char* region_name(const Pt* p, double J_T0) {
    if (p->uinf < -U_DEADBAND) {
        return (p->Vvh <= -2.0)? "descent WBS" : "vortex ring";
    }
    if (p->J > J_T0) { return "turbine"; }
    return "propeller";
}

//---- CSV ----
static void write_csv(const char* path, Pt* pts, int n) {
    FILE* f = fopen(path, "w");
    if (!f) { printf("  WARNING: could not write %s\n", path); return; }
    fprintf(f, "uinf,J,T,Q,P,CT,CP,va,Wa,a,vh,V_over_vh,mom_err,converged,admissible,ning_region\n");
    for (int i = 0; i < n; ++i) {
        Pt* p = &pts[i];
        fprintf(f, "%.4f,%.6f,%.6f,%.7f,%.4f,%.6f,%.6f,%.5f,%.5f,%.6f,%.5f,%.5f,%.8f,%d,%d,%d\n",
                p->uinf, p->J, p->T, p->Q, p->P, p->CT, p->CP, p->va, p->Wa,
                p->a, p->vh, p->Vvh, p->mom_err, p->converged, p->admissible, p->ning_region);
    }
    fclose(f);
}

static int read_csv(const char* path, Pt* pts, int max) {
    FILE* f = fopen(path, "r");
    if (!f) { return 0; }
    char line[512];
    if (!fgets(line, sizeof line, f)) { fclose(f); return 0; }   //header
    int n = 0;
    while (n < max && fgets(line, sizeof line, f)) {
        Pt* p = &pts[n];
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d",
                   &p->uinf,&p->J,&p->T,&p->Q,&p->P,&p->CT,&p->CP,&p->va,&p->Wa,
                   &p->a,&p->vh,&p->Vvh,&p->mom_err,&p->converged,&p->admissible,
                   &p->ning_region) == 16) {
            ++n;
        }
    }
    fclose(f);
    return n;
}

//---- HTML ----
//x codes: "U" freestream, "J" advance ratio, "S" state parameter V/v_h
//y codes: "T" thrust, "Q" torque, "Wa" disk axial velocity, "a" induction
//         factor, "n" normalised induced flow va/v_h, "m" momentum error %
static void emit_series(FILE* f, const char* name, Pt* p, int n,
                        const char* xf, const char* yf, double tfloor) {
    fprintf(f, "  %s: [", name);
    int first = 1;
    for (int i = 0; i < n; ++i) {
        double x, y;

        if      (strcmp(xf,"J") == 0) { x = p[i].J; }
        else if (strcmp(xf,"S") == 0) {
            //v_h is only meaningful with real thrust behind it, and V/v_h runs
            //away as T -> 0, so the state diagram is restricted to T > tfloor
            if (p[i].T <= tfloor || p[i].vh < 1e-9) { continue; }
            if (fabs(p[i].Vvh) > 4.0) { continue; }     //beyond this v_h is degenerating
            x = p[i].Vvh;
        }
        else { x = p[i].uinf; }

        if      (strcmp(yf,"T") == 0)  { y = p[i].T; }
        else if (strcmp(yf,"Q") == 0)  { y = p[i].Q; }
        else if (strcmp(yf,"Wa") == 0) { y = p[i].Wa; }
        else if (strcmp(yf,"n") == 0)  { y = p[i].va/p[i].vh; }
        else if (strcmp(yf,"a") == 0)  {
            //a = 1 - Wa/U degenerates as U -> 0; only plot where it means something
            if (fabs(p[i].uinf) < 5.0) { continue; }
            y = p[i].a;
        }
        else                           { y = p[i].mom_err*100.0; }

        if (x != x || y != y) { continue; }
        fprintf(f, "%s[%.4f,%.5f]", first? "" : ",", x, y);
        first = 0;
    }
    fprintf(f, "],\n");
}

static void write_html(const char* path, const char* tag, Pt* cur, int ncur,
                       Pt* base, int nbase, double J_T0, double tfloor) {
    FILE* f = fopen(path, "w");
    if (!f) { printf("  WARNING: could not write %s\n", path); return; }

    fputs("<!doctype html><html><head><meta charset=\"utf-8\">\n", f);
    fprintf(f, "<title>qprop flow states - %s</title>\n", tag);
    fputs("<style>\n"
          ":root{--bg:#f6f7f9;--fg:#0e1216;--mut:#868e99;--rule:#e3e7eb;--pan:#fff;--cur:#2a78d6;--base:#eb6834;--bad:#d03b3b}\n"
          "@media(prefers-color-scheme:dark){:root{--bg:#0e1013;--fg:#f1f3f6;--mut:#6f7883;--rule:#24282d;--pan:#171a1e;--cur:#3987e5;--base:#d95926;--bad:#e66767}}\n"
          "body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}\n"
          ".w{max-width:1100px;margin:0 auto;padding:32px 20px 60px}\n"
          "h1{font-size:22px;margin:0 0 4px}\n"
          "p.sub{margin:0 0 22px;color:var(--mut);font:12.5px ui-monospace,Consolas,monospace}\n"
          ".g{display:grid;grid-template-columns:1fr 1fr;gap:16px}\n"
          "@media(max-width:860px){.g{grid-template-columns:1fr}}\n"
          ".p{background:var(--pan);border:1px solid var(--rule);border-radius:3px;padding:12px 14px 6px}\n"
          ".p h2{font-size:14px;margin:0 0 2px}\n"
          ".p .n{font:11.5px ui-monospace,Consolas,monospace;color:var(--mut);margin:0 0 6px}\n"
          "svg{display:block;width:100%;height:auto}\n"
          ".lg{font:11.5px ui-monospace,Consolas,monospace;color:var(--mut);display:flex;gap:16px;margin:6px 0 10px}\n"
          ".lg i{display:inline-block;width:14px;height:2px;vertical-align:middle;margin-right:5px}\n"
          "</style></head><body><div class=\"w\">\n", f);

    fprintf(f, "<h1>qprop flow states &mdash; %s</h1>\n", tag);
    fprintf(f, "<p class=\"sub\">hqprop 5.1x3.7x3 &middot; Clark Y &middot; %.0f rpm &middot; Uinf %.0f..%.0f m/s &middot; T=0 at J=%.4f</p>\n",
            SWEEP_RPM, SWEEP_MIN, SWEEP_MAX, J_T0);
    fputs("<div class=\"lg\"><span><i style=\"background:var(--cur)\"></i>this run</span>", f);
    if (nbase) { fputs("<span><i style=\"background:var(--base)\"></i>baseline</span>", f); }
    fputs("<span><i style=\"background:var(--mut)\"></i>momentum theory</span>"
          "<span><i style=\"background:var(--bad)\"></i>inadmissible (sign(Wa) != sign(U))</span></div>\n", f);
    fputs("<div class=\"g\" id=\"g\"></div>\n<script>\nvar D={\n", f);

    emit_series(f, "T_cur",   cur, ncur, "U", "T",  tfloor);
    emit_series(f, "Q_cur",   cur, ncur, "U", "Q",  tfloor);
    emit_series(f, "n_cur",   cur, ncur, "S", "n",  tfloor);
    emit_series(f, "Wa_cur",  cur, ncur, "U", "Wa", tfloor);
    emit_series(f, "a_cur",   cur, ncur, "J", "a",  tfloor);
    emit_series(f, "m_cur",   cur, ncur, "J", "m",  tfloor);
    if (nbase) {
        emit_series(f, "T_base",  base, nbase, "U", "T",  tfloor);
        emit_series(f, "Q_base",  base, nbase, "U", "Q",  tfloor);
        emit_series(f, "n_base",  base, nbase, "S", "n",  tfloor);
        emit_series(f, "Wa_base", base, nbase, "U", "Wa", tfloor);
        emit_series(f, "a_base",  base, nbase, "J", "a",  tfloor);
        emit_series(f, "m_base",  base, nbase, "J", "m",  tfloor);
    }
    //points that fail the streamtube check, marked on the Wa panel
    fputs("  bad: [", f);
    for (int i = 0, k = 0; i < ncur; ++i) {
        if (!cur[i].admissible) { fprintf(f, "%s[%.4f,%.5f]", k++? "," : "", cur[i].uinf, cur[i].Wa); }
    }
    fputs("]\n};\n", f);

    fputs(
    "var NS='http://www.w3.org/2000/svg';\n"
    "function el(t,a){var n=document.createElementNS(NS,t);for(var k in a)n.setAttribute(k,a[k]);return n;}\n"
    "function ext(s){var lo=1e30,hi=-1e30;s.forEach(function(p){if(p[1]<lo)lo=p[1];if(p[1]>hi)hi=p[1];});return [lo,hi];}\n"
    "function chart(title,note,series,marks,xl,yl,zero){\n"
    " var W=520,H=300,ML=54,MR=16,MT=10,MB=36,PW=W-ML-MR,PH=H-MT-MB;\n"
    " var xs=[],ys=[];series.forEach(function(s){s.d.forEach(function(p){xs.push(p[0]);ys.push(p[1]);});});\n"
    " if(marks)marks.forEach(function(p){xs.push(p[0]);ys.push(p[1]);});\n"
    " var x0=Math.min.apply(null,xs),x1=Math.max.apply(null,xs);\n"
    " var y0=Math.min.apply(null,ys),y1=Math.max.apply(null,ys);\n"
    " var pad=(y1-y0)*0.08||1;y0-=pad;y1+=pad;if(zero&&y0>0)y0=0;if(zero&&y1<0)y1=0;\n"
    " var X=function(v){return ML+(v-x0)/(x1-x0)*PW;},Y=function(v){return MT+(1-(v-y0)/(y1-y0))*PH;};\n"
    " var d=document.createElement('div');d.className='p';\n"
    " d.innerHTML='<h2>'+title+'</h2><p class=\"n\">'+note+'</p>';\n"
    " var svg=el('svg',{viewBox:'0 0 '+W+' '+H});\n"
    " for(var t=0;t<=4;t++){var yv=y0+(y1-y0)*t/4;\n"
    "  svg.appendChild(el('line',{x1:ML,y1:Y(yv),x2:ML+PW,y2:Y(yv),stroke:'var(--rule)','stroke-width':1}));\n"
    "  var L=el('text',{x:ML-8,y:Y(yv)+3.5,'text-anchor':'end',fill:'var(--mut)','font-size':10,'font-family':'ui-monospace,monospace'});\n"
    "  L.textContent=yv.toFixed(Math.abs(y1-y0)<2?2:1);svg.appendChild(L);}\n"
    " for(var t=0;t<=4;t++){var xv=x0+(x1-x0)*t/4;\n"
    "  svg.appendChild(el('line',{x1:X(xv),y1:MT,x2:X(xv),y2:MT+PH,stroke:'var(--rule)','stroke-width':1}));\n"
    "  var L2=el('text',{x:X(xv),y:MT+PH+15,'text-anchor':'middle',fill:'var(--mut)','font-size':10,'font-family':'ui-monospace,monospace'});\n"
    "  L2.textContent=xv.toFixed(2);svg.appendChild(L2);}\n"
    " if(y0<0&&y1>0)svg.appendChild(el('line',{x1:ML,y1:Y(0),x2:ML+PW,y2:Y(0),stroke:'var(--mut)','stroke-width':1.25}));\n"
    " if(x0<0&&x1>0)svg.appendChild(el('line',{x1:X(0),y1:MT,x2:X(0),y2:MT+PH,stroke:'var(--mut)','stroke-width':1.25,'stroke-dasharray':'3 3'}));\n"
    " series.forEach(function(s){if(!s.d.length)return;var p='';\n"
    "  s.d.forEach(function(q,i){p+=(i?'L':'M')+X(q[0]).toFixed(2)+' '+Y(q[1]).toFixed(2);});\n"
    "  svg.appendChild(el('path',{d:p,fill:'none',stroke:s.c,'stroke-width':s.w||2,'stroke-dasharray':s.dash||'','stroke-linejoin':'round'}));});\n"
    " if(marks)marks.forEach(function(q){svg.appendChild(el('circle',{cx:X(q[0]),cy:Y(q[1]),r:2.2,fill:'var(--bad)'}));});\n"
    " var xt=el('text',{x:ML+PW/2,y:H-4,'text-anchor':'middle',fill:'var(--mut)','font-size':10,'font-family':'ui-monospace,monospace'});\n"
    " xt.textContent=xl;svg.appendChild(xt);\n"
    " var yt=el('text',{x:11,y:MT+PH/2,'text-anchor':'middle',fill:'var(--mut)','font-size':10,'font-family':'ui-monospace,monospace',transform:'rotate(-90 11 '+(MT+PH/2)+')'});\n"
    " yt.textContent=yl;svg.appendChild(yt);\n"
    " d.appendChild(svg);document.getElementById('g').appendChild(d);}\n", f);

    //split here: ISO C99 only guarantees 4095-character string literals
    fputs(
    "function S(n,c,dash,w){return {d:D[n]||[],c:c,dash:dash||'',w:w||2};}\n"
    //momentum-theory branches for the state diagram, over the data's own x-range
    "var nd=(D.n_cur||[]).concat(D.n_base||[]);\n"
    "var xlo=-4,xhi=2;\n"
    "if(nd.length){var xv=nd.map(function(p){return p[0];});\n"
    " xlo=Math.min.apply(null,xv);xhi=Math.max.apply(null,xv);}\n"
    "if(!(xhi>xlo)){xlo=-4;xhi=2;}\n"
    "D.mtProp=[];D.mtWind=[];\n"
    "for(var k=0;k<=240;k++){var v=xlo+(xhi-xlo)*k/240;\n"
    " D.mtProp.push([v,(-v+Math.sqrt(v*v+4))/2]);\n"
    " if(v<=-2)D.mtWind.push([v,-v/2-Math.sqrt(v*v/4-1)]);}\n"
    "chart('Thrust','T vs freestream velocity',[S('T_base','var(--base)','5 4'),S('T_cur','var(--cur)')],null,'Uinf (m/s)','T (N)',1);\n"
    "chart('Torque','Q vs freestream velocity',[S('Q_base','var(--base)','5 4'),S('Q_cur','var(--cur)')],null,'Uinf (m/s)','Q (N-m)',1);\n"
    "chart('Normalised induced flow','dashed grey = momentum theory; no valid branch for -2 &lt; V/v_h &lt; 0',"
    "[S('mtProp','var(--mut)','5 4',1.5),S('mtWind','var(--mut)','2 3',1.5),S('n_base','var(--base)','5 4'),S('n_cur','var(--cur)')],"
    "null,'V / v_h','va / v_h',1);\n"
    "chart('Streamtube admissibility','Wa must keep the sign of Uinf below V/v_h = -2',[S('Wa_base','var(--base)','5 4'),S('Wa_cur','var(--cur)')],D.bad,'Uinf (m/s)','Wa (m/s)',1);\n"
    "chart('Ning induction factor','a=1-Wa/U; a&gt;1 is the propeller brake region (|Uinf|&ge;5 only)',[S('a_base','var(--base)','5 4'),S('a_cur','var(--cur)')],null,'J','a',1);\n"
    "chart('Annulus-momentum error','|T_be - T_annulus| scaled by max(|T|, 5% of static T)',[S('m_base','var(--base)','5 4'),S('m_cur','var(--cur)')],null,'J','error (%)',1);\n"
    "</script></div></body></html>\n", f);
    fclose(f);
}

int main(int argc, char** argv) {
    const char* tag = (argc > 1)? argv[1] : "current";
    int strict = (getenv("STRICT") != NULL);

    Airfoil* clark_y = load_clark_y();
    if (!clark_y) { printf("TEST 8 - FAILED :( could not load polars\n"); return 1; }
    Rotor* rot = import_rotor_geometry_uiuc("../validation/hqprop_5.1x3.7x3/hqprop_5.1x3.7x3_geom.txt", clark_y, 0.12954, 3);
    if (!rot) { printf("TEST 8 - FAILED :( could not load geometry\n"); free_airfoil(clark_y); return 1; }

    double rho = 1.225, mu = 1.81e-5, a_sound = 340.0, tol = 1e-6;
    int itmax = 100;
    double Omega = SWEEP_RPM*PI/30.0;

    static Pt pts[SWEEP_N];
    int nfail_conv = 0;
    for (int i = 0; i < SWEEP_N; ++i) {
        double u = SWEEP_MIN + (SWEEP_MAX - SWEEP_MIN)*i/(double)(SWEEP_N - 1);
        solve_point(&pts[i], rot, u, Omega, rho, mu, a_sound, tol, itmax);
        if (!pts[i].converged) { ++nfail_conv; }
    }

    //thrust-reversal advance ratio, for labelling
    double J_T0 = 0.0;
    for (int i = 1; i < SWEEP_N; ++i) {
        if (pts[i-1].converged && pts[i].converged && pts[i-1].T > 0 && pts[i].T <= 0) {
            double t = pts[i-1].T/(pts[i-1].T - pts[i].T);
            J_T0 = pts[i-1].J + t*(pts[i].J - pts[i-1].J);
            break;
        }
    }

    printf("TEST 8: flow-state verification  [tag: %s%s]\n", tag, strict? ", STRICT" : "");
    printf("  rotor hqprop 5.1x3.7x3, Clark Y, %.0f rpm, Uinf %.0f..%.0f m/s, %d points\n",
           SWEEP_RPM, SWEEP_MIN, SWEEP_MAX, SWEEP_N);
    printf("  thrust reversal at J = %.4f\n\n", J_T0);

    //normalise the momentum error against the static thrust, so that points
    //near the thrust reversal do not divide a small residual by a smaller T
    double T_static = 0.0;
    for (int i = 0; i < SWEEP_N; ++i) {
        if (pts[i].converged && fabs(pts[i].uinf) < 1e-9) { T_static = fabs(pts[i].T); }
    }
    if (T_static < 1e-9) { T_static = 1.0; }
    for (int i = 0; i < SWEEP_N; ++i) {
        double scale = fabs(pts[i].T);
        if (scale < MOM_FLOOR*T_static) { scale = MOM_FLOOR*T_static; }
        pts[i].mom_err = pts[i].mom_abs/scale;
    }

    //---- per-region tally ----
    const char* regions[4] = {"descent WBS", "vortex ring", "propeller", "turbine"};
    int cnt[4] = {0,0,0,0}, adm_n[4] = {0,0,0,0}, bad_adm[4] = {0,0,0,0};
    double worst_mom[4] = {0,0,0,0};
    for (int i = 0; i < SWEEP_N; ++i) {
        if (!pts[i].converged) { continue; }
        const char* rn = region_name(&pts[i], J_T0);
        int k = 0;
        for (int j = 0; j < 4; ++j) { if (strcmp(rn, regions[j]) == 0) { k = j; break; } }
        ++cnt[k];
        if (pts[i].admissible != ADM_NA) {
            ++adm_n[k];
            if (pts[i].admissible == ADM_VIOLATED) { ++bad_adm[k]; }
        }
        if (pts[i].mom_err > worst_mom[k]) { worst_mom[k] = pts[i].mom_err; }
    }

    printf("  %-13s %6s   %-26s %s\n", "region", "points", "streamtube admissible", "momentum closure");
    printf("  ---------------------------------------------------------------------------\n");
    int total_adm = 0;
    for (int j = 0; j < 4; ++j) {
        if (!cnt[j]) { continue; }
        char adm[40];
        if (!adm_n[j]) { snprintf(adm, sizeof adm, "%-18s", "n/a"); }
        else { snprintf(adm, sizeof adm, "%3d/%-3d %-10s", adm_n[j]-bad_adm[j], adm_n[j],
                        bad_adm[j]? "<-- FAIL" : "ok"); }
        printf("  %-13s %6d   %-26s %5.1f%% worst\n", regions[j], cnt[j], adm, worst_mom[j]*100.0);
        total_adm += bad_adm[j];
    }
    printf("\n  streamtube admissibility is asserted only below V/v_h = -2, where the\n");
    printf("  windmill brake branch applies and sign(Wa) must equal sign(Uinf). In slow\n");
    printf("  descent Wa > 0 is correct, and inside the vortex ring band nothing can be\n");
    printf("  asserted. Momentum closure is reported, not asserted: Drela's circulation\n");
    printf("  closure is not expected to satisfy simple annulus momentum exactly, but a\n");
    printf("  Ning-style solver enforces it, so this column should collapse after a swap.\n");

    //---- Ning region occupancy ----
    int nreg[3] = {0,0,0};
    for (int i = 0; i < SWEEP_N; ++i) {
        if (pts[i].converged && pts[i].a == pts[i].a) { ++nreg[pts[i].ning_region]; }
    }
    printf("\n  Ning region occupancy (a = 1 - Wa/U):\n");
    printf("    momentum        a < 0.4   %4d points\n", nreg[0]);
    printf("    empirical/Buhl  0.4-1.0   %4d points\n", nreg[1]);
    printf("    propeller brake a > 1.0   %4d points\n", nreg[2]);

    //---- smoothness ----
    double max_jump = 0.0, at_J = 0.0;
    for (int i = 1; i < SWEEP_N; ++i) {
        if (!pts[i].converged || !pts[i-1].converged) { continue; }
        double d = fabs(pts[i].T - pts[i-1].T);
        if (d > max_jump) { max_jump = d; at_J = pts[i].J; }
    }
    printf("\n  smoothness: largest step in T between adjacent points = %.4f N at J = %.3f  %s\n",
           max_jump, at_J, (max_jump > TOL_SMOOTH)? "<-- FAIL" : "ok");

    //---- propeller-regime regression, via the public API ----
    printf("\n  propeller-regime regression (public qprop() entry point):\n");
    double regr_u[4] = {0.0, 10.0, 20.0, 30.0};
    int regr_bad = 0;
    for (int k = 0; k < 4; ++k) {
        RotorPerformance* perf = qprop(rot, regr_u[k], Omega, tol, itmax, rho, mu, a_sound);
        if (!perf) { printf("    Uinf %+5.1f   qprop() returned NULL  <-- FAIL\n", regr_u[k]); ++regr_bad; continue; }
        printf("    Uinf %+5.1f   T = %10.6f N   Q = %10.7f N-m   CT = %8.6f   CP = %8.6f\n",
               regr_u[k], perf->T, perf->Q, perf->CT, perf->CP);
        free_rotor_performance(perf);
    }
    printf("    (recorded in the CSV; compare against the baseline run after changing the solver)\n");

    //---- outputs ----
    char csv[256], html[256];
    snprintf(csv,  sizeof csv,  "08_flow_states_%s.csv",  tag);
    snprintf(html, sizeof html, "08_flow_states_%s.html", tag);
    write_csv(csv, pts, SWEEP_N);

    static Pt base[SWEEP_N*2];
    int nbase = 0;
    if (strcmp(tag, "before") != 0) { nbase = read_csv("08_flow_states_before.csv", base, SWEEP_N*2); }
    write_html(html, tag, pts, SWEEP_N, base, nbase, J_T0, MOM_FLOOR*T_static);

    const char* basenote = (strcmp(tag, "before") == 0)
        ? "  (this run is the baseline)"
        : (nbase? "  (baseline overlaid)" : "  (no 08_flow_states_before.csv found; run with tag 'before' first)");
    printf("\n  wrote %s\n  wrote %s%s\n", csv, html, basenote);

    //---- verdict ----
    int hard = (nfail_conv > 0 || regr_bad > 0);
    int soft = (total_adm > 0 || max_jump > TOL_SMOOTH);
    printf("\n");
    if (nfail_conv) { printf("  %d of %d sweep points failed to converge\n", nfail_conv, SWEEP_N); }
    if (total_adm)  { printf("  %d points below V/v_h = -2 violate streamtube admissibility\n", total_adm); }

    free_rotor(rot);
    free_airfoil(clark_y);

    if (hard) { printf("\nTEST 8 - FAILED :(\n"); return 1; }
    if (soft && strict) { printf("\nTEST 8 - FAILED :( (STRICT)\n"); return 1; }
    if (soft) { printf("\nTEST 8 - RECORDED (physics checks failing; set STRICT=1 to make them binding)\n"); return 0; }
    printf("\nTEST 8 - PASSED :)\n");
    return 0;
}
