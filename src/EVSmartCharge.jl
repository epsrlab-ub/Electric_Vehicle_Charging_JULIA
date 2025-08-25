module EVSmartCharge
using JuMP, Ipopt, CSV, DataFrames, LinearAlgebra, Random

export run_model, set_parameters_user, prompt_mode

# <<< paste your prompt helpers + set_parameters_user() here >>>
# <<< paste your prompt_mode() here >>>
# =========================
# Charging mode selector
# =========================
function prompt_mode()
    while true
        println("\nSelect charging type:")
        println("  1: Uncoordinated charging (no voltage constraints in scheduling)")
        println("  2: Smart charging (active power only)")
        println("  3: Smart + Reactive power")
        println("  4: Smart V2G (active discharge; no reactive)")
        println("  5: Smart V2G + Reactive power")
        print("Enter 1/2/3/4/5: ")
        s = chomp(readline())
        if !isempty(s)
            try
                m = parse(Int, s)
                if 1 <= m <= 5
                    return m
                end
            catch
                # keep looping on parse error
            end
        end
        println("Please enter 1, 2, 3, 4, or 5 (no default).")
    end
end

# ---- call it (keep this on its own line) ----
mode = prompt_mode()
@assert mode in 1:5 "Mode must be 1, 2, 3, 4, or 5"

function run_model(; data_path::String="data/", mode::Union{Int,Nothing}=nothing)
    # Load CSVs from data/
    file  = Matrix(DataFrame(CSV.File(joinpath(data_path, "InitialAdmittanceMatrix.csv"))))
    Ybus  = parse.(ComplexF64, file)
    G     = real.(Ybus); B = imag.(Ybus)

    BaseLoad  = Matrix(DataFrame(CSV.File(joinpath(data_path, "Pload_matrix.csv"))))
    BaseQLoad = Matrix(DataFrame(CSV.File(joinpath(data_path, "Qload_matrix.csv"))))
    Solar     = Matrix(DataFrame(CSV.File(joinpath(data_path, "PV_matrix.csv"))))

    # Get EV params interactively
    T, n, η, E, Δt, S_ev, SOC_initial, SOC_target, arrival_time, departure_time, ev_bus =
        set_parameters_user()
    V = n

    # Make G,B per-unit and scale loads the same way you already do
    Vpu = 12660.0; Spu = 100e6; Ypu = Spu/(Vpu^2)
    G ./= Ypu; B ./= Ypu
    BaseLoad  = (BaseLoad  .* 1000.0) ./ Spu
    BaseQLoad = (BaseQLoad .* 1000.0) ./ Spu
    Solar     = (Solar     .* 1000.0) ./ Spu
    Gen = 12660.0 / Vpu

    # If no mode passed, ask the user (no default)
    m = mode === nothing ? prompt_mode() : mode

    # <<< paste your big if/elseif block for modes 1..5 here,
    #     using the local variables above (G, B, BaseLoad, BaseQLoad, Solar, Gen, Spu, etc.) >>>

    return nothing
end

end # module
