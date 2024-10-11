# Code to merge files for different time periods
# Code for combining the time series into single nc file
# also add the lat lon values to the merged file
for ihr in range(1,25):
    datadir = f'H:/nanditha/IDF/hour{hr[ihr-1]}'
    files = sorted(glob(f'{datadir}/rcp85hotter_end*.nc'))
    combined = None
    for y in range(1,31):
        ds = xr.open_dataset(files[y-1])
        if combined is None:
            combined = ds
        else:
            combined = xr.merge([combined,ds])
        #print(y-1)
    combined.attrs = {
        "long_name": "Incremental rainfall",
        "units": "mm"
    }
    file_path = "H:/nanditha/IDF/tgw_wrf_historical_hourly_1989-09-24_01_00_00.nc"
    org_ds = xr.open_dataset(file_path)
    lat_data = org_ds.XLAT[1,:,1].values # 299 values
    lon_data = org_ds.XLONG[1,1,:].values # 424 values
    combined = combined.assign_coords(lat = lat_data, lon = lon_data)
    out_path = 'H:/nanditha/IDF/merged_files/rcp85hotter_end_hour'+str(hr[ihr-1])+'.nc'
    combined.to_netcdf(out_path)
    print(ihr-1)
    print(f'saved the file for hour {hr[ihr-1]}') 



###### lat lon in conical projections, for each latitude, longitudes spread out as we move north
import pandas as pd
import xarray as xr

# Define an empty DataFrame to store the lat/lon values
df = pd.DataFrame(columns=['lat', 'lon'])

# Load the WRF dataset
file_path = "H:/nanditha/IDF/tgw_wrf_historical_hourly_1989-09-24_01_00_00.nc"
org_ds = xr.open_dataset(file_path)

# Extract lat and lon data from the WRF file
lat_data = org_ds.XLAT[1,:,1].values  # Extract latitude data (1D array of latitudes)
lon_data = org_ds.XLONG[1,:,:].values  # Extract longitude data (2D array of longitudes)

# Initialize an empty list to collect all data
data_list = []

# Loop through each latitude and append corresponding longitude values
for ilat in range(len(lat_data)):
    # Create a dictionary for each row with the current latitude and its corresponding longitudes
    new_row = {'lat': [lat_data[ilat]] * lon_data.shape[1],  # Repeat latitude value
               'lon': lon_data[ilat, :]}  # Get the corresponding longitude values
    
    # Append this new_row dictionary to the list
    data_list.append(pd.DataFrame(new_row))
    print(f"Processed latitude index: {ilat}")

# Concatenate all the rows into a single DataFrame
df = pd.concat(data_list, ignore_index=True)

# Now 'df' contains the lat/lon grid
print(df.head())
