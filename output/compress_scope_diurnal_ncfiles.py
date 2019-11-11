#! /usr/bin/python

import sys
import numpy as np
import netCDF4
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num


# Command line arguments 

#if len(sys.argv) > 1:
#   for arg in sys.argv:
#      print 'Command line argument:', arg
#else:
#   print 'Require parameter number as command line argument... exiting'
#   sys.exit() 

# SET SOME VARIABLES

#path_to_file = '/home/563/ajn563/short/fluoro2/output/par%s/' % str(sys.argv[1])
#param = str(sys.argv[1])    # '16' 
nsplits = 20
dayint = 32   # from control file
ndays = 365   # total number of days across simulation period 
base = datetime(2010,1,1,0)


## Create datetime array for nc file

print 'Setting some datetime variables:'
print '    Start year',datetime.date(base).year
print '    Start month',datetime.date(base).month
print '    Start day',datetime.date(base).day

simulated_period = pd.date_range(datetime.date(base),periods=ndays,freq='D')

# Assign the DOY index for the 1st of every month
first_idom = np.array([0,31,59,90,120,151,181,212,243,273,304,334])

print '    hello, I hope these runs used monthly mean forcing for GPP'
print '    ...setting the datetime array to be first day of each month'
arr = np.array([base + timedelta(days=i) for i in range(ndays)])

date_arr = []

for i in range(len(first_idom)):
    date_arr.append(pd.date_range(str(arr[first_idom][i]),freq='H',periods=24))

date_arr = np.datetime_as_string(date_arr,timezone='UTC')
date_arr = date_arr.flatten().astype('datetime64')
dates_pd = pd.DatetimeIndex(date_arr)

## Loop over experiment folders and files to compress

exp_numbers = ['1000lores']
#exp_numbers_a = np.arange(1089,1093)
#exp_numbers_b = np.arange(1057,1074)
#exp_numbers_c = np.arange(1077,1089)

#exp_numbers = np.concatenate((exp_numbers_a,exp_numbers_b,exp_numbers_c))

for exp in exp_numbers:
    path = '/home/563/ajn563/short/fluoro2/output/par%s/par%s_scope_out_diurnal.nc' % (str(exp),str(exp))
    
    ## IMPORT FILE
    
    f = netCDF4.Dataset('%s' % path,'r')
    
    fgpp = f.variables['GPP'][:]
    fsif = f.variables['SIF'][:]
    ftime = f.variables['time'][:]
    ftime_units = f.variables['time'].units
    ftime_cal = f.variables['time'].calendar
    print '    Time in nc file is:',f.variables['time'].units
    vegp = f.variables['vegp'][:]
    vtype = f.variables['vtype'][:]
    gridp = f.variables['gridp'][:]
    frac = f.variables['frac'][:]
    
    f.close()
    
    
    ## Select time slices to keep
    
    ikeep = np.squeeze(np.where(fgpp[0,:] > -990))
    
    fsif_out = fsif[:,ikeep]
    fgpp_out = fgpp[:,ikeep]
    
    ## Write sliced data and datetime to netcdf file
    
    nvp = len(vegp)
    
    path_outputncf = '/home/563/ajn563/short/fluoro2/output/par%s/par%s_scope_out.nc' % (str(exp),str(exp))

    dataset = netCDF4.Dataset('%s' % path_outputncf, 'w', format='NETCDF4_CLASSIC')
    
    # Create Dimensions
    dataset.createDimension('time',None)
    dataset.createDimension('vegp',nvp)
    
    # Create Variables
    nc_times = dataset.createVariable('time', np.float64, ('time',))
    nc_vegp = dataset.createVariable('vegp', np.int32, ('vegp',))
    nc_gridp = dataset.createVariable('gridp', np.int32, ('vegp',))
    nc_frac = dataset.createVariable('frac', np.float64, ('vegp',))
    nc_vtype = dataset.createVariable('vtype', np.int32, ('vegp',))
    nc_sif = dataset.createVariable('SIF', np.float64, ('vegp','time'))
    nc_gpp = dataset.createVariable('GPP', np.float64, ('vegp','time'))
    
    # Create some Global Attributes
    dataset.description = 'BETHY-SCOPE model simulation of SIF and GPP'
    dataset.source = 'Model runs performed on NCI (Raijin) by Alex Norton'
    
    # Create some Variable Attributes
    nc_times.units = 'hours since 0001-01-01 00:00:00.0'
    nc_times.calendar = 'gregorian'
    nc_vegp.long_name = 'Model vegetation point'
    nc_gridp.long_name = 'Model grid point'
    nc_frac.long_name = 'Grid point fractional coverage of vegp'
    nc_vtype.long_name = 'Vegetation type'
    nc_sif.long_name = 'SCOPE Solar-Induced Fluorescence at 755 nm'
    nc_sif.units = 'W m-2 um-1 sr-1'
    nc_gpp.long_name = 'SCOPE Gross Primary Production'
    nc_gpp.units = 'umol C m-2 s-1'
    
    # Fill variables with data
    nc_times[:] = date2num(dates_pd.astype(datetime), units=ftime_units, calendar=ftime_cal)
    nc_vegp[:] = vegp
    nc_gridp[:] = gridp
    nc_frac[:] = frac
    nc_vtype[:] = vtype
    nc_sif[:,:] = fsif_out
    nc_gpp[:,:] = fgpp_out
    
    dataset.close()




