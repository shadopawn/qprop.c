#-------------------------------------------------------------------------------
#   Validation Script: Groupner 6x3 Propeller Aerodynamic Analysis
#
#   This script validates the accuracy of qprop.c by comparing its predictions
#   with values returned by the original QPROP (v1.22)
#   The present validation case uses the Groupner 6x3 propeller, the default
#   case proposed by QPROP.
#
#   How to run:
#   python3 graupner6x3_validation_Uinf=5.py
#
#   Author: Andrea Pavan
#   License: MIT
#-------------------------------------------------------------------------------
import ctypes;
import math;
import matplotlib.pyplot as plt;
import os;
import sys;
sys.path.insert(0, "./qprop-portable/");
import qprop;

def main():
    #define airfoil polars using the analytical model of the original QPROP
    #the coefficients are matching those written in cam6x3.def
    airfoil_analytic = qprop.analytic_polar_curves(
        0.50, 5.8, -0.3, 1.2,           #CL0, CL_a, CLmin, CLmax
        0.028, 0.050, 0.020, 0.5,       #CD0, CD2u, CD2l, CLCD0
        70000.0, -0.7                   #REref, REexp
    )

    #read propeller geometry data from the original QPROP output
    original_output = []
    with open(os.path.join("original_qprop1.22_data_Uinf=5","cam6x3_qprop1.22_output.txt"), "r") as fileio:
        lines = fileio.readlines()[24:]     #skip the first 24 lines
        for line in lines:
            original_output.append([float(x) for x in line.split()])
    #the variable original_output is a List of List
    #original_output[i][1]: radial distance of the element centers (m)
    #original_output[i][2]: chord of each element (m)
    #original_output[i][3]: sweep angle of each element (deg)

    #create propeller geometry
    #NOTE that the output file radius is from the element center, not from the section
    sections = []
    dr_elem = original_output[1][0] - original_output[0][0]     #size of first element
    r_elem = original_output[0][0]                              #location of first element center
    r = r_elem - dr_elem/2                                      #location of root station
    dc_elem = original_output[1][1] - original_output[0][1]     #chord variation in the first element
    c_elem = original_output[0][1]
    c = c_elem - dc_elem
    dβ_elem = original_output[1][2] - original_output[0][2]     #twist variation in the first element
    β_elem = original_output[0][2]
    β = β_elem - dβ_elem
    sections.append(qprop.create_section(
        c,
        qprop.deg2rad(β),
        r,
        airfoil_analytic
    ))
    for i in range(1, len(original_output)):
        r_elem = original_output[i][0]
        c_elem = original_output[i][1]
        β_elem = original_output[i][2]
        if i >= 1:
            dr_elem = r_elem - original_output[i-1][0]
            dc_elem = c_elem - original_output[i-1][1]
            dβ_elem = β_elem - original_output[i-1][2]
        r = r_elem + dr_elem/2
        c = c_elem + dc_elem/2
        β = β_elem + dβ_elem/2
        sections.append(qprop.create_section(
            c,
            qprop.deg2rad(β),
            r,
            airfoil_analytic
        ))
    D = 2 * sections[-1].r
    B = 2
    graupner6x3 = qprop.create_rotor(D, B, len(sections), sections)

    #run qprop.c
    Uinf = 5.0;                         #freestream velocity (m/s)
    Omega = 14020*math.pi/30;           #rotor speed (rad/s)
    qpropc_results = qprop.qprop(graupner6x3, Uinf, Omega, 1e-6, 200)
    nelems = qpropc_results.nelems
    for i in range(nelems):
        if abs(qpropc_results.residuals[i]) > 1e-6:
            print("ERROR while running qprop: convergence not reached in one or more elements")
            break
    print("qprop.c results:")
    print("  Thrust: ", round(qpropc_results.T, 5), " N")
    print("  Torque: ", round(qpropc_results.Q, 5), " N-m")

    #compare with original QPROP results
    r_original = [row[0] for row in original_output]
    dr_original = [0.0] * len(r_original)
    dr_original[0] = r_original[1] - r_original[0]
    for i in range(1, len(r_original)-1):
        dr_original[i] = 0.5 * (r_original[i+1] - r_original[i-1])
    dr_original[-1] = r_original[-1] - r_original[-2]
    Wa_original = [row[9] for row in original_output]
    Wt_original = [Wa * row[0] / (row[11] * (D/2))                  for Wa, row in zip(Wa_original, original_output)]
    W_original = [math.sqrt(Wa**2 + Wt**2)                          for Wa, Wt in zip(Wa_original, Wt_original)]
    phi_original = [math.atan(Wa/Wt)                                for Wa, Wt in zip(Wa_original, Wt_original)]
    Cl_original = [row[3] for row in original_output]
    Cd_original = [row[4] for row in original_output]
    Cn_original = [Cl * math.cos(phi) - Cd * math.sin(phi)          for Cl, Cd, phi in zip(Cl_original, Cd_original, phi_original)]
    Ct_original = [Cl * math.sin(phi) + Cd * math.cos(phi)          for Cl, Cd, phi in zip(Cl_original, Cd_original, phi_original)]
    dTdr_original = [0.5 * 1.225 * W**2 * Cn * row[1]               for W, Cn, row in zip(W_original, Cn_original, original_output)]
    dQdr_original = [0.5 * 1.225 * W**2 * Ct * row[1] * row[0]      for W, Ct, row in zip(W_original, Ct_original, original_output)]
    print("QPROP v1.22 results:")
    print(f"  Thrust:{round(2 * sum(dTdr_original[i] * dr_original[i] for i in range(len(dr_original))), 5)}N")
    print(f"  Torque:{round(2 * sum(dQdr_original[i] * dr_original[i] for i in range(len(dr_original))), 5)}N-m")

    #compare thrust distributions
    plt1 = plt.figure(figsize=(6,4), dpi=100)       #600x400px
    plt.plot(
        [qpropc_results.r[i] / (D/2)        for i in range(nelems)],
        [qpropc_results.dTdr[i]             for i in range(nelems)],
        label = "qprop.c",
        linewidth = 2
    )
    plt.scatter(
        [r_original[i] / (D/2)              for i in range(nelems)],
        [dTdr_original[i]                   for i in range(nelems)],
        label = "QPROP v1.22",
        marker = "D",
        s = 20,
        color = "darkorange"
    )
    plt.title("Graupner 6x3 Thrust (Uinf=5.0m/s)")
    plt.xlabel("Blade radius r/R")
    plt.ylabel("Thrust distribution dT/dr (N/m)")
    plt.grid(True, which="both", linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.show()

    #compare torque distributions
    plt2 = plt.figure(figsize=(6,4), dpi=100)       #600x400px
    plt.plot(
        [qpropc_results.r[i] / (D/2)        for i in range(nelems)],
        [qpropc_results.dQdr[i]             for i in range(nelems)],
        label = "qprop.c",
        linewidth = 2
    )
    plt.scatter(
        [r_original[i] / (D/2)              for i in range(nelems)],
        [dQdr_original[i]                   for i in range(nelems)],
        label = "QPROP v1.22",
        marker = "D",
        s = 20,
        color = "darkorange"
    )
    plt.title("Graupner 6x3 Torque (Uinf=5.0m/s)")
    plt.xlabel("Blade radius r/R")
    plt.ylabel("Torque distribution dQ/dr (N-m/m)")
    plt.grid(True, which="both", linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
