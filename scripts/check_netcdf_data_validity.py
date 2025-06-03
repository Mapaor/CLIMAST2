import numpy as np
from netCDF4 import Dataset

# --------- CONFIGURATION ---------
variable = 'tas'  # 'tas' for temperature or 'pr' for precipitation
membre = 'em' # 'r1', 'r2'... 'r32' per triar un membre, o 'em' per ensemble mean
ensemble_file = f"./member_data/{variable}_Amon_IPSL-CM6A-LR_historical_{membre}i1p1f1_gr_185001-201412.nc"

# --------- LOAD ENSEMBLE MEAN DATA ---------
try:
    with Dataset(ensemble_file, mode='r') as nc:
        ensemble_mean = nc.variables[variable][:]
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        time = nc.variables['time'][:]
except FileNotFoundError:
    print(f"File {ensemble_file} not found. Check the path and filename.")
    exit()

# --------- CHECK DATA ---------
print(f"Checking data for variable: {variable}")
print(f"Shape of ensemble mean data: {ensemble_mean.shape}")
print(f"Longitude range: {lon.min()} to {lon.max()}")
print(f"Latitude range: {lat.min()} to {lat.max()}")
print(f"Time range: {time.min()} to {time.max()}")

# Compute basic statistics
data_min = np.min(ensemble_mean)
data_max = np.max(ensemble_mean)
data_mean = np.mean(ensemble_mean)

print(f"Minimum value: {data_min}")
print(f"Maximum value: {data_max}")
print(f"Mean value: {data_mean}")

# Check if values are within expected ranges
if variable == 'tas':  # Temperature in Celsius
    if not (0 <= data_min <= 373 and 0 <= data_max <= 373):
        print("Warning: Temperature values are outside the expected range (-100 to 100°C).")
elif variable == 'pr':  # Precipitation in mm/day
    if not (0 <= data_min <= 500 and 0 <= data_max <= 500):
        print("Warning: Precipitation values are outside the expected range (0 to 500 mm/day).")

print("Data check complete.")

### RESULTATS ESPERATS (ordre's de magnitud):
#
# --- Temperatura ---
# Checking data for variable: tas
# Shape of ensemble mean data: (1980, 143, 144)
# Longitude range: 0.0 to 357.5
# Latitude range: -90.0 to 90.0
# Time range: 15.5 to 60249.5
# Minimum value: 198.10545349121094
# Maximum value: 312.1200256347656
# Mean value: 276.7853088378906
# Data check complete.
#
# --- Precipitació ---
# Checking data for variable: pr
# Shape of ensemble mean data: (1980, 143, 144)
# Longitude range: 0.0 to 357.5
# Latitude range: -90.0 to 90.0
# Time range: 15.5 to 60249.5
# Minimum value: 0.0
# Maximum value: 0.0011421270901337266
# Mean value: 2.821780799422413e-05
# Data check complete.