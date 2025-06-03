# ---------------------------------------------------
#        Aquest script genera un fitxer NetCDF 
#          amb l'ensemble mean a partir dels
#          fitxers dels 32 membres del model.
#        Nota: pot tardar fins a uns 10 minuts
# ---------------------------------------------------

import numpy as np
from netCDF4 import Dataset
import os

# --------- CONFIGURATION ---------
variable = 'pr'  # 'tas' for temperature or 'pr' for precipitation
num_members = 32
source_data = './member_data/'
output_file = f"./member_data/{variable}_Amon_IPSL-CM6A-LR_historical_emi1p1f1_gr_185001-201412.nc"

# --------- LOAD DATA ---------
fields = []
for i in range(1, num_members + 1):
    file_name = f"{source_data}{variable}_Amon_IPSL-CM6A-LR_historical_r{i}i1p1f1_gr_185001-201412.nc"
    try:
        with Dataset(file_name, mode='r') as fh:
            field = fh.variables[variable][:]  # Load data without unit conversion
            fields.append(field)
            if i == 1:  # Save dimensions and coordinates from the first file
                lon = fh.variables['lon'][:]
                lat = fh.variables['lat'][:]
                time = fh.variables['time'][:]
    except FileNotFoundError:
        print(f"File {file_name} not found. Skipping.")

# --------- COMPUTE ENSEMBLE MEAN ---------
fields = np.array(fields)
ensemble_mean = np.mean(fields, axis=0)

# --------- SAVE TO NETCDF ---------
with Dataset(output_file, mode='w', format='NETCDF4') as nc_out:
    # Create dimensions
    nc_out.createDimension('time', len(time))
    nc_out.createDimension('lat', len(lat))
    nc_out.createDimension('lon', len(lon))

    # Create variables
    times = nc_out.createVariable('time', 'f4', ('time',))
    lats = nc_out.createVariable('lat', 'f4', ('lat',))
    lons = nc_out.createVariable('lon', 'f4', ('lon',))
    ensemble_var = nc_out.createVariable(variable, 'f4', ('time', 'lat', 'lon'), zlib=True)

    # Assign data to variables
    times[:] = time
    lats[:] = lat
    lons[:] = lon
    ensemble_var[:, :, :] = ensemble_mean

    # Add attributes
    times.units = 'days since 1850-01-01'
    lats.units = 'degrees_north'
    lons.units = 'degrees_east'
    ensemble_var.units = 'kg/m²/s' if variable == 'pr' else 'Kelvin'
    ensemble_var.long_name = f"Ensemble mean of {variable}"

print(f"Ensemble mean saved to {output_file}")