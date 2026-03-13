#-------------------------------------------------------------------------------
#   Validation Script: Groupner 6x3 Propeller Aerodynamic Analysis
#
#   This script validates the accuracy of qprop.c by comparing its predictions
#   with values returned by the original QPROP (v1.22)
#   The present validation case uses the Groupner 6x3 propeller, the default
#   case proposed by QPROP.
#
#   How to run:
#   julia graupner6x3_validation_hovering.jl
#
#   Author: Andrea Pavan
#   License: MIT
#-------------------------------------------------------------------------------
using DelimitedFiles;
using Plots;
include("qprop-portable/qprop.jl");
import .QProp;

function main()
    #define airfoil polars using the analytical model of the original QPROP
    #the coefficients are matching those written in cam6x3.def
    airfoil_analytic = QProp.analytic_polar_curves(
        0.50, 5.8, -0.3, 1.2,           #CL0, CL_a, CLmin, CLmax
        0.028, 0.050, 0.020, 0.5,       #CD0, CD2u, CD2l, CLCD0
        70000.0, -0.7                   #REref, REexp
    );

    #read propeller performance data from the original QPROP output
    original_output = readdlm(
        joinpath(@__DIR__, "original_qprop1.22_data_hovering", "cam6x3_qprop1.22_output.txt"),
        Float64,
        header = false,
        skipstart = 24
    );
    #the variable original_output is a 25×12 Matrix{Float64}:
    #original_output[:,1]: radial distance of the element centers (m)
    #original_output[:,2]: chord of each element (m)
    #original_output[:,3]: sweep angle of each element (deg)

    #create propeller geometry
    #NOTE that the output file radius is from the element center, not from the section
    sections = Vector{QProp.Section}(undef, 0);
    dr_elem = original_output[2,1] - original_output[1,1];      #size of first element
    r_elem = original_output[1,1];                              #location of first element center
    r = r_elem - dr_elem/2;                                     #location of root station
    dc_elem = original_output[2,2] - original_output[1,2];      #chord variation in the first element
    c_elem = original_output[1,2];
    c = c_elem - dc_elem;
    dβ_elem = original_output[2,3] - original_output[1,3];      #twist variation in the first element
    β_elem = original_output[1,3];
    β = β_elem - dβ_elem;
    push!(sections, QProp.Section(
        c,
        deg2rad(β),
        r,
        airfoil_analytic
    ));
    for i in 1:lastindex(original_output,1)
        #create ending section for the i-th element
        r_elem = original_output[i,1];
        c_elem = original_output[i,2];
        β_elem = original_output[i,3];
        if i>=2
            dr_elem = r_elem - original_output[i-1,1];
            dc_elem = c_elem - original_output[i-1,2];
            dβ_elem = β_elem - original_output[i-1,3];
        end
        r = r_elem + dr_elem/2;
        c = c_elem + dc_elem/2;
        β = β_elem + dβ_elem/2;
        push!(sections, QProp.Section(
            c,
            deg2rad(β),
            r,
            airfoil_analytic
        ));
    end
    D = 2 * sections[end].r;
    B = 2;
    graupner6x3 = QProp.Rotor(D, B, length(sections), sections);
    
    #ALTERNATIVE: create propeller geometry from the original QPROP input file
    #=sections = Vector{QProp.Section}(undef, 0);
    open(joinpath(@__DIR__, "original_qprop1.22_data_hovering", "cam6x3.def"), "r") do fileio
        line_counter = 0;
        for line in eachline(fileio)
            line = strip(line);
            if isempty(line)
                continue;
            end
            if !startswith(line,'#') && !startswith(line,'!')
                line_counter += 1;
            end
            if line_counter < 9
                continue;
            end
            
            #read propeller geometry
            line_fields = split(line);
            r = parse(Float64, line_fields[1]);
            c = parse(Float64, line_fields[2]);
            β = parse(Float64, line_fields[3]);
            push!(sections, QProp.Section(
                0.0254 * c,
                deg2rad(β),
                0.0254 * r,
                airfoil_analytic
            ));
        end
    end
    D = 2 * sections[end].r;
    B = 2;
    graupner6x3 = QProp.Rotor(D, B, length(sections), sections);
    nsections = 50;
    graupner6x3 = QProp.refine_rotor_sections(graupner6x3, nsections);=#

    #run qprop.c
    Uinf = 0.01;                #freestream velocity (m/s)
    Ω = 14020*pi/30;            #rotor speed (rad/s)
    qpropc_results = QProp.qprop(graupner6x3, Uinf, Ω, 1e-6, 200);
    if any(abs.(qpropc_results.residuals) .> 1e-6)
        error("ERROR while running qprop: convergence not reached in one or more elements");
    end
    println("qprop.c results:");
    println("  Thrust: ", round(qpropc_results.T, digits=5), " N");
    println("  Torque: ", round(qpropc_results.Q, digits=5), " N-m");

    #compare with original QPROP results
    r_original = original_output[:,1];
    dr_original = similar(r_original);
    dr_original[1] = r_original[2] - r_original[1];
    dr_original[2:end-1] = 0.5*(r_original[3:end] - r_original[1:end-2]);
    dr_original[end,1] = r_original[end] - r_original[end-1];
    c_original = original_output[:,2];
    Wa_original = original_output[:,10];
    Wt_original = Wa_original .* r_original ./ (original_output[:,12] * (D/2));
    W_original = sqrt.(Wa_original.^2 + Wt_original.^2);
    phi_original = atan.(Wa_original./Wt_original);
    #vt_original = Wa_original .* tan.(deg2rad.(original_output[:,11]));
    Cl_original = original_output[:,4];
    Cd_original = original_output[:,5];
    Cn_original = Cl_original.*cos.(phi_original) - Cd_original.*sin.(phi_original);
    Ct_original = Cl_original.*sin.(phi_original) + Cd_original.*cos.(phi_original);
    dTdr_original = 0.5 * 1.225 * W_original.^2 .* Cn_original .* c_original;
    dQdr_original = 0.5 * 1.225 * W_original.^2 .* Ct_original .* c_original .* r_original;
    println("QPROP v1.22 results:");
    println("  Thrust: ", round(B*sum(dTdr_original.*dr_original), digits=5), " N");
    println("  Torque: ", round(B*sum(dQdr_original.*dr_original), digits=5), " N-m");

    #compare thrust distributions
    plt1 = plot(qpropc_results.r/(D/2), qpropc_results.dTdr, label="qprop.c", linewidth=2,
        title = "Graupner 6x3 Thrust (Hovering)",
        xlabel = "Blade radius r/R",
        ylabel = "Thrust distribution dT/dr (N/m)",
        minorgrid = true
    );
    scatter!(plt1, r_original/(D/2), dTdr_original, label="QPROP v1.22", markershape=:diamond, markersize=4);
    display(plt1);

    #compare torque distributions
    plt2 = plot(qpropc_results.r/(D/2), qpropc_results.dQdr, label="qprop.c", linewidth=2,
        title = "Graupner 6x3 Torque (Hovering)",
        xlabel = "Blade radius r/R",
        ylabel = "Torque distribution dQ/dr (N-m/m)",
        minorgrid = true
    );
    scatter!(plt2, r_original/(D/2), dQdr_original, label="QPROP v1.22", markershape=:diamond, markersize=4);
    display(plt2);
end

main();
