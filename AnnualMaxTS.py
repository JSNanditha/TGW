'''
Array Job

Code to generate nc file for each year with  maximum precipitation for each water year for durations 1 hour to 24 hours
'''
## loading libraries

import xarray as xr
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from glob import glob
import os

print("Libraries loaded")

# reading files
datadir = '/scratch/Users/njayadevan/rcp85cooler_2060_2099/hourly/'
files = sorted(glob('{}*.nc'.format(datadir)))


print('working')
# Create a date range from 1989-10-01 to 2018-10-01
time_coords = pd.date_range(start='1989-10-01', end='2019-10-01', freq='Y')
#time_coords = pd.date_range(start='2029-10-01', end='2059-10-01', freq='Y')
#time_coords = pd.date_range(start='2069-10-01', end='2099-10-01', freq='Y')
# Convert to numpy array if needed
time_coords = time_coords.to_numpy()
#years = time_coords.astype('datetime64[Y]').astype(int) + 1970
years = time_coords.astype('datetime64[Y]').astype(int)+1970
print(years)
# Get the job index (1-based) from the environment variable
y = int(os.getenv('SGE_TASK_ID'))  # Replace SGE_TASK_ID with the appropriate variable for your scheduler if needed
#
#y=1
# indexing for the first date, i.e. October 1st 1989
ind1 = 506 + (y-1)*52
out = None
    # reading 52 files for each water year
for ifile in range(ind1,ind1+53): #559

	 # opening the first entry in the file
    ds = xr.open_dataset(files[ifile])
    sel_var = ds[['RAINC', 'RAINSH','RAINNC','Time']]
    rainc = sel_var['RAINC']
    rainsh = sel_var['RAINSH']
    rainnc = sel_var['RAINNC']
   
    rain = rainc+rainnc+rainsh
    rain.attrs = {
        "long_name": "Accumulated rainfall",
        "units": "mm"
    }


     # appending all the files for each year (i.e. 52 + 1 previous week file)   
    if out is None:
        out = rain
    else:
        out = xr.concat([out,rain],dim = 'Time')
    
    print(files[ifile])
lind = len(out)-1
    ## estimating a rolling sum for different durations

for itime in range(1,25):
	# taking sum from October 1st to Sep 30 of each water year
	#24*7 = 168
	# 168 hours in Sep 24 file
	# after taking difference, the data range reduces by 1
	# we need to remove the first 167 hours which is the iter of 166, we exclude it and count from 167
  # beacuse last week from Sep need to be removed

    # rolling sum from duration 1 to 24 after subtracting the previous day from the current day, then estimate the maximum value for each water year
    roll_sum = out.diff(dim = 'Time').rolling(Time=itime, center=True).sum()[167:lind].max(dim = 'Time')
    
# adding the year to each nc file and changing the name of the variable to 'pr'
    
    roll_sum = roll_sum.expand_dims({'Time': [years[y-1]]}).rename('pr')
    roll_sum.attrs = {
     "long_name": "Incremental rainfall",
        "units": "mm"
    }
# creating a directory for each hour and storing all the variables inside
    #os.makedirs(path+'hour'+str(itime),exist_ok = True)
    output_path = path+'hour'+str(itime)+'/rcp85cooler_end_ts_'+str(years[y-1])+'.nc'
    roll_sum.to_netcdf(output_path)
    
print(y)
print("Saved file for year "+str(years[y-1])+" at "+str(output_path)+ 'for all durations')
