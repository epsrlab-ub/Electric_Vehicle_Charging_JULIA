module EVSmartCharge
using JuMP, Ipopt, CSV, DataFrames, LinearAlgebra, Random

export run_model, set_parameters_user, prompt_mode

# -------------------------
# User parameter utilities
# -------------------------

"""
    _prompt_default(msg, default)

Utility used by `set_parameters_user` to query the user for a value while
providing a sensible default. When running in a non-interactive environment the
default is returned immediately, avoiding blocking on `readline`.
"""
function _prompt_default(msg::String, default)
    if !isatty(stdin)
        return default
    end
    print("$msg [$default]: ")
    s = try
        chomp(readline())
    catch
        ""
    end
    return isempty(s) ? default : parse(typeof(default), s)
end

"""
    set_parameters_user()

Collect basic EV and simulation parameters from the user. Reasonable defaults
are provided so that the function can run unattended during automated tests.
"""
function set_parameters_user()
    T            = _prompt_default("Number of time periods", 24)
    n            = _prompt_default("Number of buses", 33)
    η            = _prompt_default("Charging/discharging efficiency", 0.95)
    E            = _prompt_default("Battery capacity (kWh)", 60.0)
    Δt           = _prompt_default("Time step (h)", 1.0)
    S_ev         = _prompt_default("EV power rating (kW)", 7.0)
    SOC_initial  = _prompt_default("Initial SOC (0-1)", 0.2)
    SOC_target   = _prompt_default("Target SOC (0-1)", 0.9)
    arrival_time = _prompt_default("Arrival time step", 1)
    departure_time = _prompt_default("Departure time step", T)
    ev_bus       = _prompt_default("Bus hosting the EV", 2)
    return T, n, η, E, Δt, S_ev, SOC_initial, SOC_target, arrival_time, departure_time, ev_bus
end

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
    # ------------------------------------------------------------------
    # Helper for data loading (handles files with a trailing `1`)
    # ------------------------------------------------------------------
    function _load_matrix(fname)
        path = joinpath(data_path, fname)
        if !isfile(path)
            # fall back to a version with "1" inserted before extension
            base, ext = splitext(fname)
            path = joinpath(data_path, base * "1" * ext)
        end
        Matrix(DataFrame(CSV.File(path)))
    end

    # Load CSVs from data/
    file  = _load_matrix("InitialAdmittanceMatrix.csv")
    Ybus  = parse.(ComplexF64, file)
    G     = real.(Ybus); B = imag.(Ybus)

    BaseLoad  = _load_matrix("Pload_matrix.csv")
    BaseQLoad = _load_matrix("Qload_matrix.csv")
    Solar     = _load_matrix("PV_matrix.csv")

    # Get EV params interactively
    T, n, η, E, Δt, S_ev, SOC_initial, SOC_target, arrival_time, departure_time, ev_bus =
        set_parameters_user()

    # Dimensions implied by data (override user if inconsistent)
    n = size(BaseLoad, 1)
    T = size(BaseLoad, 2)

    # Per-unit conversion for network parameters
    Vpu = 12_660.0
    Spu = 100e6
    Ypu = Spu / (Vpu^2)
    G ./= Ypu
    B ./= Ypu
    BaseLoad  = (BaseLoad  .* 1000.0) ./ Spu
    BaseQLoad = (BaseQLoad .* 1000.0) ./ Spu
    Solar     = (Solar     .* 1000.0) ./ Spu

    # Pre-compute base voltage drop from loads and solar
    const_drop = zeros(n, T)
    for t in 1:T, v in 1:n
        for u in 1:n
            const_drop[v, t] += G[v, u] * (BaseLoad[u, t] - Solar[u, t]) +
                                B[v, u] * BaseQLoad[u, t]
        end
    end

    # If no mode passed, ask the user (no default)
    m = mode === nothing ? prompt_mode() : mode

    # -----------------------
    # Build optimization model
    # -----------------------
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    # Common variables
    @variable(model, 0 <= soc[1:T] <= 1)
    @variable(model, V[1:n, 1:T])

    # Active power variable bounds depend on mode
    if m in (1,2,3)
        @variable(model, 0 <= p[1:T] <= S_ev)
    else
        @variable(model, -S_ev <= p[1:T] <= S_ev)
    end

    # Reactive power variable for modes with Q support
    if m in (3,5)
        @variable(model, -S_ev <= q[1:T] <= S_ev)
    end

    # Fix power to zero when EV is not connected
    for t in 1:T
        if t < arrival_time || t >= departure_time
            @constraint(model, p[t] == 0)
            if m in (3,5)
                @constraint(model, q[t] == 0)
            end
        end
    end

    # Uncoordinated charging: force max charging whenever connected
    if m == 1
        for t in arrival_time:min(departure_time-1, T)
            @constraint(model, p[t] == S_ev)
        end
    end

    # SOC dynamics
    @constraint(model, soc[1] == SOC_initial)
    for t in 1:T-1
        expr = η * p[t] * Δt / E
        @constraint(model, soc[t+1] == soc[t] + expr)
    end
    @constraint(model, soc[T] >= SOC_target)

    # Voltage equations and limits
    for t in 1:T, v in 1:n
        if m in (3,5)
            @constraint(model, V[v, t] == 1 - const_drop[v, t] -
                                      G[v, ev_bus] * p[t] -
                                      B[v, ev_bus] * q[t])
        else
            @constraint(model, V[v, t] == 1 - const_drop[v, t] -
                                      G[v, ev_bus] * p[t])
        end
        if m != 1
            @constraint(model, 0.95 <= V[v, t] <= 1.05)
        end
    end

    # Objective: simple quadratic charging cost + reactive effort
    obj = @expression(model, sum(p[t]^2 for t in 1:T))
    if m in (3,5)
        obj += @expression(model, sum(q[t]^2 for t in 1:T))
    end
    @objective(model, Min, obj)

    optimize!(model)

    result = Dict(
        :p => value.(p),
        :soc => value.(soc),
        :voltage => value.(V),
    )
    if m in (3,5)
        result[:q] = value.(q)
    end
    return result
end

end # module
