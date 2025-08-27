#!/usr/bin/env julia
# Entry point script for running the EVSmartCharge model

# Make the source directory available on the load path and import the module
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using EVSmartCharge

function main()
    try
        result = run_model()
        println("Model run completed.")
        if result !== nothing
            println("Result: ", result)
        end
    catch err
        println("Model run failed: ", err)
    end
end

main()
