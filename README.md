# Electric_Vehicle_Charging_JULIA
This project implements a distribution network + electric vehicle (EV) charging optimization framework using Julia. The model is built with JuMP and solved with Ipopt, allowing multiple charging strategies to be simulated and compared on a per-bus basis.
## Contents
- **Uncoordinated Charging**
- **Smart Charging**
- **Smart Charging + Reactive Power Support**
- **Vehicle-to-grid**
- **Vehicle-to-grid + Reactive Power Support**
## Contributions

This work is the result of collaborative efforts:

| Name                 | Role                      | Contact                |
|----------------------|---------------------------|------------------------|
| **Jordan Ngamaleu-Kemta** | Developed the set of codes | 📧 Jgngamal@buffalo.edu@buffalo.edu |
| **Adedoyin Inaolaji** | Supervised the work        | 📧 ainaolaj@buffalo.edu |
## Features

- **IEEE 33-bus distribution feeder** with hourly load and solar generation data  
- **EV fleet scheduling** with arrival/departure randomness and SOC targets  
- **Five charging modes**, from uncoordinated charging to V2G with reactive power  
- **Per-unit system** (100 MVA, 12.66 kV base) for grid consistency  
- **Power-flow solved** via current-injection model  
- **Voltage limits enforced** (0.95–1.05 p.u.)  
- **Plots of bus voltages and EV SOC trajectories**
