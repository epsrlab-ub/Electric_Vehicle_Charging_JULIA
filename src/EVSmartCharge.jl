module EVSmartCharge
using JuMP, Ipopt, CSV, DataFrames, LinearAlgebra, Random

export run_model, set_parameters_user, prompt_mode
# ------------------------------------------------------------
# Helper functions for interactive input
# ------------------------------------------------------------

"""Safely read an integer from `readline()`.

Continues prompting until the user provides a valid integer. Optional
`min` and `max` bounds can be supplied to restrict the range of valid
values.
"""
function read_int(prompt; min::Union{Nothing,Int}=nothing, max::Union{Nothing,Int}=nothing)
    while true
        print(prompt)
        s = chomp(readline())
        try
            v = parse(Int, s)
            if (min !== nothing && v < min) || (max !== nothing && v > max)
                println("Value must be between $(min === nothing ? "-∞" : string(min)) and $(max === nothing ? "∞" : string(max)).")
            else
                return v
            end
        catch
            println("Please enter a valid integer.")
        end
    end
end

"""Safely read a floating point number from `readline()`.

Behavior mirrors [`read_int`] but parses `Float64` values.
"""
function read_float(prompt; min::Union{Nothing,Float64}=nothing, max::Union{Nothing,Float64}=nothing)
    while true
        print(prompt)
        s = chomp(readline())
        try
            v = parse(Float64, s)
            if (min !== nothing && v < min) || (max !== nothing && v > max)
                println("Value must be between $(min === nothing ? "-∞" : string(min)) and $(max === nothing ? "∞" : string(max)).")
            else
                return v
            end
        catch
            println("Please enter a valid number.")
        end
    end
end

"""Read a vector of numbers from `readline()`.

`count` specifies how many entries are expected. The `parse_fn` argument
is applied elementwise to convert the strings (e.g. `parse.(Int, ...)`).
"""
function read_vector(prompt, count, parse_fn)
    while true
        print(prompt)
        parts = split(chomp(readline()))
        if length(parts) != count
            println("Please enter $count values separated by spaces.")
            continue
        end
        try
            return [parse_fn(p) for p in parts]
        catch
            println("Invalid input; try again.")
        end
    end
end

read_vector_int(prompt, count) = read_vector(prompt, count, x -> parse(Int, x))
read_vector_float(prompt, count) = read_vector(prompt, count, x -> parse(Float64, x))

# ------------------------------------------------------------
# Collect EV parameters from the user
# ------------------------------------------------------------
function set_parameters_user()
    println("\nEnter simulation parameters (press Enter after each value).")
    T  = read_int("Number of time periods T: "; min=1)
    Δt = read_float("Length of each period Δt (hours): "; min=0.0)
    η  = read_float("Charging efficiency η (0-1): "; min=0.0, max=1.0)

    # Buses that host EVs
    nbus = read_int("Number of buses with EVs: "; min=0)
    buses = nbus == 0 ? Int[] : read_vector_int("Enter the bus indices (space-separated): ", nbus)

    E = Float64[]            # Battery capacity (kWh) per EV
    S_ev = Float64[]         # Max charging rate (kW) per EV
    SOC_initial = Float64[]  # Arrival state of charge (0-1)
    SOC_target  = Float64[]  # Target SOC at departure (default 1.0)
    arrival_time = Int[]     # Arrival time step
    departure_time = Int[]   # Departure time step
    ev_bus = Int[]           # Bus index for each EV

    for bus in buses
        nev = read_int("\nNumber of EVs at bus $bus: "; min=1)
        for ev in 1:nev
            push!(ev_bus, bus)
            push!(E, read_float("  EV $ev battery capacity (kWh): "; min=0.0))
            push!(S_ev, read_float("  EV $ev max charging rate (kW): "; min=0.0))
            push!(SOC_initial, read_float("  EV $ev arrival SOC (0-1): "; min=0.0, max=1.0))
            push!(SOC_target, 1.0)  # assume full charge desired
            at = read_int("  EV $ev arrival time (1-$T): "; min=1, max=T)
            dt = read_int("  EV $ev departure time ($at-$T): "; min=at, max=T)
            push!(arrival_time, at)
            push!(departure_time, dt)
        end
    end

    n = length(ev_bus)
    return (T, n, η, E, Δt, S_ev, SOC_initial, SOC_target, arrival_time, departure_time, ev_bus)
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

"""
Run a simple EV charging optimization.

The model schedules charging power for each EV so that the desired
state-of-charge is reached by the specified departure time while
minimizing total charging power. Only very basic constraints are
considered, but this function serves as a concrete example of how the
collected user inputs can be fed into an optimization model.
"""
function run_model()
    T, n, η, E, Δt, S_ev, SOC_initial, SOC_target, arrival_time, departure_time, ev_bus = set_parameters_user()
    mode = prompt_mode()

    model = Model(Ipopt.Optimizer)

    @variable(model, 0 <= SOC[1:T, 1:n] <= 1)
    @variable(model, 0 <= P[1:T-1, 1:n])

    for v in 1:n
        @constraint(model, SOC[1, v] == SOC_initial[v])
        for t in 1:(T-1)
            if arrival_time[v] <= t < departure_time[v]
                @constraint(model, P[t, v] <= S_ev[v])
                @constraint(model, SOC[t+1, v] == SOC[t, v] + (η * P[t, v] * Δt) / E[v])
            else
                @constraint(model, P[t, v] == 0)
                @constraint(model, SOC[t+1, v] == SOC[t, v])
            end
        end
        @constraint(model, SOC[departure_time[v], v] >= SOC_target[v])
    end

    @objective(model, Min, sum(P))
    optimize!(model)

    return (mode=mode, SOC=value.(SOC), P=value.(P))
end

end # module EVSmartCharge
