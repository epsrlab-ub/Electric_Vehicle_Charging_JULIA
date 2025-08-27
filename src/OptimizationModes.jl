module OptimizationModes

using JuMP, Ipopt

# Simple integer prompt helper used by the mode selector
function prompt_int(msg, default)
    print(msg, ": ")
    line = readline()
    return try
        parse(Int, line)
    catch
        default
    end
end

# =========================
# Charging mode selector
# =========================
function prompt_mode()
    println("\nSelect charging type:")
    println("  1: Uncoordinated charging (no voltage constraints in scheduling)")
    println("  2: Smart charging (active power only)")
    println("  3: Smart + Reactive power")
    println("  4: Smart V2G (active discharge; no reactive)")
    println("  5: Smart V2G + Reactive power")
    m = prompt_int("Enter 1/2/3/4/5", 1)
    @assert m in (1,2,3,4,5) "Mode must be 1, 2, 3, 4, or 5"
    return m
end

# ---------------------------------------------------------------------------
# Top-level driver for the optimization modes. The function expects all
# necessary data (network matrices, EV parameters, etc.) to be defined in the
# calling scope. This mirrors the original script used for experimentation and
# can be adapted into the main application as needed.
# ---------------------------------------------------------------------------
function run_optimization()
    mode = prompt_mode()

    if mode == 1
        # -----------------------------------------
        # Mode 1: Uncoordinated per-EV scheduling (no grid constraints),
        #         then run a power-flow with fixed EV load.
        # -----------------------------------------
        n_ev = V
        SOC_results        = zeros(T, n_ev)
        P_charging_results = zeros(T-1, n_ev)   # P defined on 1..T-1

        for ev in 1:n_ev
            UNCOORD = Model(Ipopt.Optimizer)

            @variable(UNCOORD, 0 <= SOC_ev[1:T] <= 1)
            @variable(UNCOORD, 0 <= P_ev_uc[1:T-1] <= S_ev)

            @constraint(UNCOORD, SOC_ev[1] == SOC_initial[ev])

            for t in 1:(T-1)
                if (arrival_time[ev] <= t) && (t < departure_time[ev])
                    @constraint(UNCOORD, SOC_ev[t+1] == SOC_ev[t] + (η * P_ev_uc[t] * Δt) / E)
                else
                    @constraint(UNCOORD, P_ev_uc[t] == 0.0)
                    @constraint(UNCOORD, SOC_ev[t+1] == SOC_ev[t])
                end
            end

            @objective(UNCOORD, Min, sum((SOC_ev[t] - SOC_target[ev])^2 for t in 1:T))
            optimize!(UNCOORD)

            SOC_results[:, ev]        .= value.(SOC_ev)
            P_charging_results[:, ev] .= value.(P_ev_uc)
        end

        # Aggregate EV kW -> per-bus, per-time (per-unit for PF)
        Ev_Load_uc = zeros(nbus, T)  # p.u.
        for t in 1:(T-1), v in 1:V
            Ev_Load_uc[ev_bus[v], t] += (1000.0/Spu) * P_charging_results[t, v]
        end
        # t = T stays 0.0 by definition

        # ---------- Power Flow (current-injection model) ----------
        PF = Model(Ipopt.Optimizer)

        @variable(PF, Pgen[1:nbus, 1:T])
        @variable(PF, Qgen[1:nbus, 1:T])
        @variable(PF, V_Real[1:nbus, 1:T], start = Gen)
        @variable(PF, V_Imag[1:nbus, 1:T], start = 0.0)
        @variable(PF, I_Real[1:nbus, 1:T])
        @variable(PF, I_Imag[1:nbus, 1:T])

        # I = Y * V (rectangular)
        for t in 1:T
            @constraint(PF, [i=1:nbus], I_Real[i,t] == sum(G[i,j]*V_Real[j,t] - B[i,j]*V_Imag[j,t] for j in 1:nbus))
            @constraint(PF, [i=1:nbus], I_Imag[i,t] == sum(B[i,j]*V_Real[j,t] + G[i,j]*V_Imag[j,t] for j in 1:nbus))

            # Power balance with fixed EV load (per-unit)
            @NLconstraint(PF, [i=1:nbus],
                Pgen[i,t] - BaseLoad[i,t] + Solar[i,t] - Ev_Load_uc[i,t]== V_Real[i,t]*I_Real[i,t] + V_Imag[i,t]*I_Imag[i,t])
            @NLconstraint(PF, [i=1:nbus],Qgen[i,t] - BaseQLoad[i,t]== V_Imag[i,t]*I_Real[i,t] - V_Real[i,t]*I_Imag[i,t])
        end

        # Slack & generator settings
        @constraint(PF, [t=1:T], V_Real[1,t] == Gen)
        @constraint(PF, [t=1:T], V_Imag[1,t] == 0.0)
        @constraint(PF, [i=2:nbus, t=1:T], Pgen[i,t] == 0.0)
        @constraint(PF, [i=2:nbus, t=1:T], Qgen[i,t] == 0.0)

        # Minimize active power losses
        Ploss = @expression(PF, sum(Pgen[i,t] - BaseLoad[i,t] + Solar[i,t] - Ev_Load_uc[i,t] for i in 1:nbus, t in 1:T))
        @objective(PF, Min, Ploss)

        optimize!(PF)

        V_Real_uc = value.(V_Real); V_Imag_uc = value.(V_Imag)
    elseif mode == 2
        # -----------------------------------------
        # Mode 2: Smart (active power only)
        # -----------------------------------------
        SMARTCHARGE = Model(Ipopt.Optimizer)

        @variable(SMARTCHARGE, Pgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, Qgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, V_Real[1:nbus, 1:T], start = Gen)
        @variable(SMARTCHARGE, V_Imag[1:nbus, 1:T], start = 0.0)
        @variable(SMARTCHARGE, I_Real[1:nbus, 1:T])
        @variable(SMARTCHARGE, I_Imag[1:nbus, 1:T])

        @variable(SMARTCHARGE, 0 <= SOC[1:T, 1:V] <= 1)
        @variable(SMARTCHARGE, 0 <= P_ev[1:(T-1), 1:V])
        @constraint(SMARTCHARGE, [t=1:T-1,v=1:V], P_ev[t,v] <= S_ev)

        for v in 1:V
            @constraint(SMARTCHARGE, SOC[1, v] == SOC_initial[v])
            for t in 1:(T-1)
                if (arrival_time[v] <= t) && (t < departure_time[v])
                    @constraint(SMARTCHARGE, SOC[t+1, v] == SOC[t, v] + (η * P_ev[t, v] * Δt) / E)
                else
                    @constraint(SMARTCHARGE, P_ev[t, v] == 0.0)
                    @constraint(SMARTCHARGE, SOC[t+1, v] == SOC[t, v])
                end
            end
        end

        @expression(SMARTCHARGE, Ev_Load[i=1:nbus, t=1:T],t <= T-1 ? (1000/Spu) * sum(P_ev[t, v] for v in 1:V if ev_bus[v] == i) : 0.0)
        @expression(SMARTCHARGE, Q_ev_bus_pu_expr[i=1:nbus, t=1:T], 0.0)

        for t in 1:T
            for i in 1:nbus
                @constraint(SMARTCHARGE, I_Real[i, t] == sum(G[i,j]*V_Real[j,t] - B[i,j]*V_Imag[j,t] for j in 1:nbus))
                @constraint(SMARTCHARGE, I_Imag[i, t] == sum(B[i,j]*V_Real[j,t] + G[i,j]*V_Imag[j,t] for j in 1:nbus))
            end
            for i in 1:nbus
                @NLconstraint(SMARTCHARGE, Pgen[i,t] - BaseLoad[i,t] + Solar[i,t] - Ev_Load[i,t] == V_Real[i,t]*I_Real[i,t] + V_Imag[i,t]*I_Imag[i,t])
                @NLconstraint(SMARTCHARGE, Qgen[i,t] - BaseQLoad[i,t] - Q_ev_bus_pu_expr[i,t] == V_Imag[i,t]*I_Real[i,t] - V_Real[i,t]*I_Imag[i,t])
            end
        end

        @constraint(SMARTCHARGE, [t=1:T], V_Real[1,t] == Gen)
        @constraint(SMARTCHARGE, [t=1:T], V_Imag[1,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Pgen[i,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Qgen[i,t] == 0.0)
        @NLconstraint(SMARTCHARGE, [i=1:nbus, t=1:T], 0.95^2 <= V_Real[i,t]^2 + V_Imag[i,t]^2 <= 1.05^2)

        @objective(SMARTCHARGE, Min,sum( (SOC[t, v] - SOC_target[v])^2 for t in 1:T, v in 1:V ))
        optimize!(SMARTCHARGE)

    elseif mode == 3
        # -----------------------------------------
        # Mode 3: Smart + Reactive power
        # -----------------------------------------
        P_Max = S_ev
        S_app = 1.1 * P_Max
        Q_Max = sqrt(S_app^2 - P_Max^2)

        SMARTCHARGE = Model(Ipopt.Optimizer)

        @variable(SMARTCHARGE, Pgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, Qgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, V_Real[1:nbus, 1:T], start = Gen)
        @variable(SMARTCHARGE, V_Imag[1:nbus, 1:T], start = 0.0)
        @variable(SMARTCHARGE, I_Real[1:nbus, 1:T])
        @variable(SMARTCHARGE, I_Imag[1:nbus, 1:T])

        @variable(SMARTCHARGE, 0 <= SOC[1:T, 1:V] <= 1)
        @variable(SMARTCHARGE, 0 <= P_ev[1:(T-1), 1:V] <= P_Max)
        @variable(SMARTCHARGE, -Q_Max <= Q_ev[1:(T-1), 1:V] <= Q_Max)

        for v in 1:V
            @constraint(SMARTCHARGE, SOC[1, v] == SOC_initial[v])
            for t in 1:(T-1)
                if (arrival_time[v] <= t) && (t < departure_time[v])
                    @constraint(SMARTCHARGE, SOC[t+1, v] == SOC[t, v] + (η * P_ev[t, v] * Δt) / E)
                    @NLconstraint(SMARTCHARGE, P_ev[t, v]^2 + Q_ev[t, v]^2 <= S_app^2)
                else
                    @constraint(SMARTCHARGE, P_ev[t, v] == 0.0)
                    @constraint(SMARTCHARGE, Q_ev[t, v] == 0.0)
                    @constraint(SMARTCHARGE, SOC[t+1, v] == SOC[t, v])
                end
            end
        end

        @expression(SMARTCHARGE, Ev_Load[i=1:nbus, t=1:T],t <= T-1 ? (1000/Spu) * sum(P_ev[t, v] for v in 1:V if ev_bus[v] == i) : 0.0)
        @expression(SMARTCHARGE, Q_ev_bus_pu_expr[i=1:nbus, t=1:T],t <= T-1 ? (1000/Spu) * sum(Q_ev[t, v] for v in 1:V if ev_bus[v] == i) : 0.0)

        for t in 1:T
            for i in 1:nbus
                @constraint(SMARTCHARGE, I_Real[i, t] == sum(G[i,j]*V_Real[j,t] - B[i,j]*V_Imag[j,t] for j in 1:nbus))
                @constraint(SMARTCHARGE, I_Imag[i, t] == sum(B[i,j]*V_Real[j,t] + G[i,j]*V_Imag[j,t] for j in 1:nbus))
            end
            for i in 1:nbus
                @NLconstraint(SMARTCHARGE, Pgen[i,t] - BaseLoad[i,t] + Solar[i,t] - Ev_Load[i,t] == V_Real[i,t]*I_Real[i,t] + V_Imag[i,t]*I_Imag[i,t])
                @NLconstraint(SMARTCHARGE, Qgen[i,t] - BaseQLoad[i,t] - Q_ev_bus_pu_expr[i,t] == V_Imag[i,t]*I_Real[i,t] - V_Real[i,t]*I_Imag[i,t])
            end
        end

        @constraint(SMARTCHARGE, [t=1:T], V_Real[1,t] == Gen)
        @constraint(SMARTCHARGE, [t=1:T], V_Imag[1,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Pgen[i,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Qgen[i,t] == 0.0)
        @NLconstraint(SMARTCHARGE, [i=1:nbus, t=1:T], 0.95^2 <= V_Real[i,t]^2 + V_Imag[i,t]^2 <= 1.05^2)

        @objective(SMARTCHARGE, Min,sum( (SOC[t, v] - SOC_target[v])^2 for t in 1:T, v in 1:V ))
        optimize!(SMARTCHARGE)

    elseif mode == 4
        # -----------------------------------------
        # Mode 4: Smart V2G (active discharge; no reactive from EVs)
        # -----------------------------------------
        SMARTCHARGE = Model(Ipopt.Optimizer)

        @variable(SMARTCHARGE, Pgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, Qgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, V_Real[1:nbus, 1:T], start = Gen)
        @variable(SMARTCHARGE, V_Imag[1:nbus, 1:T], start = 0.0)
        @variable(SMARTCHARGE, I_Real[1:nbus, 1:T])
        @variable(SMARTCHARGE, I_Imag[1:nbus, 1:T])

        @variable(SMARTCHARGE, 0 <= SOC[1:T, 1:V] <= 1)
        @variable(SMARTCHARGE, 0 <= P_char[1:(T-1), 1:V])
        @variable(SMARTCHARGE, 0 <= P_dis[1:(T-1), 1:V])
        @constraint(SMARTCHARGE, [t=1:T-1,v=1:V], P_char[t,v] <= S_ev)
        @constraint(SMARTCHARGE, [t=1:T-1,v=1:V], P_dis[t,v] <= S_ev)

        P_net = @expression(SMARTCHARGE, [t=1:T-1, v=1:V], P_char[t,v] - P_dis[t,v])

        for v in 1:V
            @constraint(SMARTCHARGE, SOC[1, v] == SOC_initial[v])
        end

        for v in 1:V
            for t in 1:(T-1)
                if (arrival_time[v] <= t) && (t < departure_time[v])
                    @constraint(SMARTCHARGE,SOC[t+1, v] == SOC[t, v] +(η * P_char[t, v] * Δt)/E - ((1/η) * P_dis[t, v]  * Δt)/E )
                else
                    @constraint(SMARTCHARGE, P_char[t, v] == 0.0)
                    @constraint(SMARTCHARGE, P_dis[t, v] == 0.0)
                    @constraint(SMARTCHARGE, SOC[t+1, v] == SOC[t, v])
                end
            end
        end

        @expression(SMARTCHARGE, Ev_Load[i=1:nbus, t=1:T],t <= T-1 ? (1000/Spu) * sum(P_net[t, v] for v in 1:V if ev_bus[v] == i) : 0.0)

        for t in 1:T
            for i in 1:nbus
                @constraint(SMARTCHARGE,
                    I_Real[i, t] == sum(G[i,j]*V_Real[j,t] - B[i,j]*V_Imag[j,t] for j in 1:nbus))
                @constraint(SMARTCHARGE,
                    I_Imag[i, t] == sum(B[i,j]*V_Real[j,t] + G[i,j]*V_Imag[j,t] for j in 1:nbus))
            end
            for i in 1:nbus
                @NLconstraint(SMARTCHARGE,
                    Pgen[i,t] - BaseLoad[i,t] + Solar[i,t] - Ev_Load[i,t] == V_Real[i,t]*I_Real[i,t] + V_Imag[i,t]*I_Imag[i,t])
                @NLconstraint(SMARTCHARGE,
                    Qgen[i,t] - BaseQLoad[i,t] == V_Imag[i,t]*I_Real[i,t] - V_Real[i,t]*I_Imag[i,t])
            end
        end

        @constraint(SMARTCHARGE, [t=1:T], V_Real[1,t] == Gen)
        @constraint(SMARTCHARGE, [t=1:T], V_Imag[1,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Pgen[i,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Qgen[i,t] == 0.0)
        @NLconstraint(SMARTCHARGE, [i=1:nbus, t=1:T], 0.95^2 <= V_Real[i,t]^2 + V_Imag[i,t]^2 <= 1.05^2)

        @objective(SMARTCHARGE, Min, sum((SOC[t, v] - SOC_target[v])^2 for t in 1:T, v in 1:V ))

        optimize!(SMARTCHARGE)

    elseif mode == 5
        # -----------------------------------------
        # Mode 5: Smart V2G + Reactive power
        # -----------------------------------------
        P_Max = S_ev
        S_app = 1.1 * P_Max
        Q_Max = sqrt(S_app^2 - P_Max^2)

        SMARTCHARGE = Model(Ipopt.Optimizer)

        @variable(SMARTCHARGE, Pgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, Qgen[1:nbus, 1:T])
        @variable(SMARTCHARGE, V_Real[1:nbus, 1:T], start = Gen)
        @variable(SMARTCHARGE, V_Imag[1:nbus, 1:T], start = 0.0)
        @variable(SMARTCHARGE, I_Real[1:nbus, 1:T])
        @variable(SMARTCHARGE, I_Imag[1:nbus, 1:T])

        @variable(SMARTCHARGE, 0 <= SOC[1:T, 1:V] <= 1)
        @variable(SMARTCHARGE, 0 <= P_char[1:(T-1), 1:V] <= P_Max)
        @variable(SMARTCHARGE, 0 <= P_dis[1:(T-1), 1:V] <= P_Max)
        @variable(SMARTCHARGE, -Q_Max <= Q_ev[1:(T-1), 1:V] <= Q_Max)

        P_net = @expression(SMARTCHARGE, [t=1:T-1, v=1:V], P_char[t,v] - P_dis[t,v])

        for v in 1:V
            @constraint(SMARTCHARGE, SOC[1, v] == SOC_initial[v])
        end

        for v in 1:V
            for t in 1:(T-1)
                if (arrival_time[v] <= t) && (t < departure_time[v])
                    @constraint(SMARTCHARGE,SOC[t+1, v] == SOC[t, v] +(η     * P_char[t, v] * Δt)/E - ((1/η) * P_dis[t, v]  * Δt)/E )
                    @NLconstraint(SMARTCHARGE, P_net[t, v]^2 + Q_ev[t, v]^2 <= S_app^2)
                else
                    @constraint(SMARTCHARGE, P_char[t, v] == 0.0)
                    @constraint(SMARTCHARGE, P_dis[t, v] == 0.0)
                    @constraint(SMARTCHARGE, Q_ev[t, v]  == 0.0)
                    @constraint(SMARTCHARGE, SOC[t+1, v] == SOC[t, v])
                end
            end
        end

        @expression(SMARTCHARGE, Ev_Load[i=1:nbus, t=1:T],
            t <= T-1 ? (1000/Spu) * sum(P_net[t, v] for v in 1:V if ev_bus[v] == i) : 0.0)
        @expression(SMARTCHARGE, Ev_Load_rtc[i=1:nbus, t=1:T],t <= T-1 ? (1000/Spu) * sum(Q_ev[t, v]  for v in 1:V if ev_bus[v] == i) : 0.0)

        for t in 1:T
            for i in 1:nbus
                @constraint(SMARTCHARGE,
                    I_Real[i, t] == sum(G[i,j]*V_Real[j,t] - B[i,j]*V_Imag[j,t] for j in 1:nbus))
                @constraint(SMARTCHARGE,
                    I_Imag[i, t] == sum(B[i,j]*V_Real[j,t] + G[i,j]*V_Imag[j,t] for j in 1:nbus))
            end
            for i in 1:nbus
                @NLconstraint(SMARTCHARGE,
                    Pgen[i, t] - BaseLoad[i, t] + Solar[i, t] - Ev_Load[i, t] == V_Real[i,t]*I_Real[i,t] + V_Imag[i,t]*I_Imag[i,t])
                @NLconstraint(SMARTCHARGE,
                    Qgen[i, t] - BaseQLoad[i, t] - Ev_Load_rtc[i, t] == V_Imag[i,t]*I_Real[i,t] - V_Real[i,t]*I_Imag[i,t])
            end
        end

        @constraint(SMARTCHARGE, [t=1:T], V_Real[1,t] == Gen)
        @constraint(SMARTCHARGE, [t=1:T], V_Imag[1,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Pgen[i,t] == 0.0)
        @constraint(SMARTCHARGE, [i=2:nbus, t=1:T], Qgen[i,t] == 0.0)
        @NLconstraint(SMARTCHARGE, [i=1:nbus, t=1:T], 0.95^2 <= V_Real[i,t]^2 + V_Imag[i,t]^2 <= 1.05^2)

        @objective(SMARTCHARGE, Min, sum((SOC[t, v] - SOC_target[v])^2 for t in 1:T, v in 1:V ))
        optimize!(SMARTCHARGE)
    end
end

end # module OptimizationModes

