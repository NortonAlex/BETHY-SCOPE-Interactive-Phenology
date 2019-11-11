#! /usr/bin/python

import sys
import numpy as np
import netCDF4
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num

# Import my own functions
import sys
sys.path.append('/home/563/ajn563/short/pyFunctions/')
from my_pyFunctions import *

## Command line arguments
#  first argument = spatial resolution of model run (i.e. lores or hires)
#  second argument = number of blocks/splits

# Command line arguments 

if len(sys.argv) > 1:
   for arg in sys.argv:
      print 'Command line argument:', arg
else:
   print 'Require parameter number as command line argument... exiting'
   sys.exit()

## SET SOME VARIABLES

# Choose the variables to stitch together onto the lat-lon grid
variables = ['rpasm','rfluo','rgppfluo','rnpp','rnep']

exp = str(sys.argv[1])
resolution = str(sys.argv[2])
nblocks = int(sys.argv[3]) 

ndays = 365
tsteps = 12

path_to_output = '/home/563/ajn563/short/fluoro2/output/%s/' % exp
base = datetime(2015,1,1,0)    # starting datetime for simulated period
path_to_vpblockfile = '/home/563/ajn563/short/fluoro2/control_blockruns/block_vplist.txt'

print '   Output path:', path_to_output


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

# Only select one time per month as it's a monthly run
tinds = np.arange(0,12*24,24)
dates_pd = dates_pd[tinds]

## IMPORT GRID-VEG DATA

path_to_gridveg = '/home/563/ajn563/short/fluoro2/control_bethy/'

if resolution == 'lores':
    print '   low-resolution grid scale specified'
    f = netCDF4.Dataset('%sbethy_loresgrid3veg_faparratio.nc' % path_to_gridveg,'r')    # lores
    gridnum = f.variables['gridnum'][::-1,:]
    vtype = f.variables['type'][:,::-1,:]
    frac = f.variables['frac'][:,::-1,:]
    lats = f.variables['lat'][::-1]
    lons = f.variables['lon'][:]
    f.close()
    latstep = np.abs(lats[1] - lats[0])
    lonstep = np.abs(lons[1] - lons[0])
    ngp = 170
    nvp = 506
    # Per veg-point scale info
    f = netCDF4.Dataset('/home/563/ajn563/short/fluoro2/output/lores_scope_gridinfo.nc','r')
    gridp_vp = f.variables['gridp'][:]
    f.close()
elif resolution == 'hires':
    print '   high-resolution grid scale specified'
    f = netCDF4.Dataset('%sbethy_grid3veg_v2.nc' % path_to_gridveg,'r')        # hires
    gridnum = f.variables['gridnum'][::-1,:]
    vtype = f.variables['type'][:,::-1,:]
    frac = f.variables['frac'][:,::-1,:]
    lats = f.variables['lat'][::-1]
    lons = f.variables['lon'][:]
    f.close()
    latstep = np.abs(lats[1] - lats[0])
    lonstep = np.abs(lons[1] - lons[0])
    ngp = 3462
    nvp = 8776
    # Per veg-point scale info
    #f = netCDF4.Dataset('/home/563/ajn563/short/fluoro2/output/hires_scope_gridinfo.nc','r')
    #gridp_vp = f.variables['gridp'][:].data
    #f.close()  
    f = '/home/563/ajn563/short/fluoro2/output/hires_scope_gridinfo.txt'
    gridp_vp = np.loadtxt(f,dtype=int)
else:
    print '   No or incorrect model grid resolution specified. Must specify as 2nd command line argument on call... exiting'
    sys.exit()

## Create Land-Ocean grid mask for BETHY

frac_nan = np.ma.getmask(np.isnan(frac[0,:,:]))
print 'shape:', np.shape(frac_nan)
tlen = 12
frac_stacked = np.array([frac_nan for _ in range(tlen)])
print 'shape:', np.shape(frac_stacked)
oceanmask = frac_stacked
#print 'Total number of grid points:', np.size(oceanmask[0,:,:])
#print '             Ocean grid points:', np.shape(np.where(oceanmask[0,:,:] == True))
#print '              Land grid points (with vegetation):', np.shape(np.where(oceanmask[0,:,:] == False))
#print ''
#print 'Ocean Mask:'
#print '   True = ocean'
#print '  False = land (vegetation points)'

xlon, ylat = np.meshgrid(lons, lats)

igp_landpoints = np.where(oceanmask == False)
igp_landpoints_geog = np.where(oceanmask[0,:,:] == False)
print 'shape of mask:', np.shape(igp_landpoints)

lats_gps = ylat[igp_landpoints_geog]
lons_gps = xlon[igp_landpoints_geog]


# Initiate output array and dictionary
var_output = np.zeros((tsteps,ngp))
outputs_dict = {}


## Get indexes of what grid points were run per block 
igp_dict = {}

fp = open('%s' % path_to_vpblockfile,'r')

ilines = np.arange(1,nblocks*2,2)

blockcount = 0
for i, line in enumerate(fp):
    if (i % 2 == 1):
        blockcount+=1
        #print 'Even line in blockvpfile, specifying the veg-points for this block (%02i)' % blockcount
        #print '  VP from:',line[0:3],'to',line[-5:]
        line_vps = line.split()
        vps = np.array(line_vps,dtype=int)
        #print '  VP from:',vps[0], 'to', vps[-1]
        vps_inds = vps - 1
        
        i0 = vps_inds[0]
        i1 = vps_inds[-1]
        gps = np.unique(gridp_vp[i0:i1+1])    # grid-points run for this block
        gps_inds = gps - np.int(1)    # indexes of grid-points run for this block
        #print '  GP from:',gps[0],'to',gps[-1]
        igp_dict['%02i' % blockcount] = gps_inds
        
fp.close()


## Import variables of interest
## - Loop over each variable,
## - Get relevant grid-points from each block run,
## - Stitch together these blocks and,
## - Map the variable back onto lat-lon grid.


for var in variables:
    print 'Variable:',var
    blockcount = 0
    for iblock in range(nblocks):
        block = iblock+1
        
        ## Get land grid point indexes for this block
        gps_inds = igp_dict['%02i' % block]
        
        ## Import variable file
        #print '   File: %sjob%02i%s.nc' % (path_to_output,blockcount,var)
        f = netCDF4.Dataset('%sjob%02i/%s.nc' % (path_to_output,block,var),'r')
        var_data = f.variables['%s' % var][:,::-1,:]
        ## roll lon-axis 180 degrees
        var_data = np.roll(var_data,len(lons)/2,axis=2)
        f.close()
        
        ## If it's SIF or GPP, we normalize the summed rates (the 24 time-steps are summed per day in the code)
        if (var == 'rfluo') | (var == 'rgppfluo'):
            var_data = var_data/24.    # gives you the daily mean rate of gpp or sif (e.g. umol C m-2 s-1)
        
        ## Pull out land points data
        for i in range(tsteps):
            var_data_tslice = var_data[i,:,:]
            var_data_tgp = var_data_tslice[igp_landpoints_geog]    # all land grid points for this time-slice
            var_output[i,gps_inds] = var_data_tgp[gps_inds]
    ## Map field back onto Lat-Lon grid
    field_gridded = np.zeros(np.shape(oceanmask))
    #print np.shape(field_gridded)
    for t in range(12):
        maplats, maplons, data_map, n_points = maptrack(lats_gps,lons_gps,var_output[t,:],maplats=lats[::-1],maplons=lons,latres=latstep,lonres=lonstep)
        #print '  data_mapped shape:',np.shape(data_map)
        field_gridded[t,:,:] = data_map
    ## Add this field to the output dictionary
    outputs_dict['%s' % var] = field_gridded

fp.close()


## Write variables to netCDF file
## - one netcdf file per variable

for var in variables:
    
    dataset = netCDF4.Dataset('%s%s_stitched.nc' % (path_to_output,var), 'w', format='NETCDF4_CLASSIC')
    
    # Create Dimensions
    dataset.createDimension('time',None)
    dataset.createDimension('latitude',len(lats))
    dataset.createDimension('longitude',len(lons))
    
    # Create Variables
    nc_times = dataset.createVariable('time', np.float64, ('time',))
    nc_lats = dataset.createVariable('latitude', np.float64, ('latitude',))
    nc_lons = dataset.createVariable('longitude', np.float64, ('longitude',))
    nc_var = dataset.createVariable('%s' % var, np.float64, ('time','latitude','longitude'))
    
    # Set lat-lon details
    nc_lats.axis = "Y"
    nc_lons.axis = "X"
    nc_lats.units="degrees_north"
    nc_lons.units="degrees_east"
    nc_lats.standard_name = "latitude"
    nc_lons.standard_name="longitude"
    
    # Create some Global Attributes
    dataset.description = 'BETHY-SCOPE model simulation'
    dataset.source = 'Model runs performed on NCI (Raijin) by Alex Norton'
    
    # Create some Variable Attributes
    nc_times.units = 'hours since 0001-01-01 00:00:00.0'
    nc_times.calendar = 'gregorian'
    #nc_var.long_name = 'BETHY-SCOPE %s' % var
    if var == 'rpasm':
        nc_var.long_name = 'BETHY Plant Available Soil Moisture'
        nc_var.units = 'mm'
    elif var == 'rfluo':
        nc_var.long_name = 'BETHY-SCOPE Solar Induced Fluorescence (daily mean)'
        nc_var.units = 'W m-2 um-1 sr-1'
    elif var == 'rgppfluo':
        nc_var.long_name = 'BETHY-SCOPE Gross Primary Production (daily mean)'
        nc_var.units = 'umol m-2 s-1'
    elif var == 'rnpp':
        nc_var.long_name = 'BETHY-SCOPE Net Primary Production (unsure on units)'
        nc_var.units = ''
    elif var == 'rnep':
        nc_var.long_name = 'BETHY-SCOPE Net Ecosystem Production (unsure on units)'
        nc_var.units = ''
    
    # Fill variables with data
    nc_times[:] = date2num(dates_pd.astype(datetime), units=nc_times.units, calendar=nc_times.calendar)
    nc_lats[:] = maplats
    nc_lons[:] = maplons
    nc_var[:,:,:] = outputs_dict['%s' % var]
    
    dataset.close()

