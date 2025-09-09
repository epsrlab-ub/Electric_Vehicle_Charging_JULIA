#!/usr/bin/env julia
# Entry point script for running the EVSmartCharge model

# Activate the project environment and ensure dependencies are available
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Make the source directory available on the load path and import the module
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using EVSmartChargedcdcdc

function main()
    try
        result = run_model()
        println("Model run completed.")
        println("Selected mode: ", result.mode)
        println("Charging power schedule (kW):\n", result.P)
        println("State of charge trajectory:\n", result.SOC)
    catch err
        println("Model run failed: ", err)
    end
end

main()
