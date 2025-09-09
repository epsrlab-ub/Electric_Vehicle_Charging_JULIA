# Convert Power Data consisting of real and imaginary load to matrices of 24 hour time window
using DelimitedFiles, XLSX

# Load the power data
loaddata = readdlm("powerdata_PV.txt")
Pload = loaddata[:, 2]
Qload = loaddata[:, 3]
PV = loaddata[:, 4]

# Load the normalized curves from the Excel file
xlsx_data = XLSX.readxlsx("Normalized_curves.xlsx")
LoadTimes = xlsx_data[1][:, 1]
LoadTimes = transpose(LoadTimes)
PVTimes = xlsx_data[2][:, 1]
PVTimes = transpose(PVTimes)

# Creating the 33x24 Matrices
Pload_matrix = Pload * LoadTimes
Qload_matrix = Qload * LoadTimes
PV_matrix = PV * PVTimes

# Testing
# Pload_matrix = Pload * ones(1,24)
# Qload_matrix = Qload * ones(1,24)
# PV_matrix = PV * ones(1,24)

# Save the matrices to CSV files
writedlm("Pload_matrix.csv", Pload_matrix, ',')
writedlm("Qload_matrix.csv", Qload_matrix, ',')
writedlm("PV_matrix.csv", PV_matrix, ',')
