#! /usr/bin/python

import sys
import numpy as np
import netCDF4
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num


# Command line arguments 

if len(sys.argv) > 1:
   for arg in sys.argv:
      print 'Command line argument:', arg
else:
   print 'Require parameter number as command line argument... exiting'
   sys.exit() 

# SET SOME VARIABLES

path_to_output = '/home/563/ajn563/short/fluoro2/output/par%s/' % str(sys.argv[1])
param = str(sys.argv[1])    # '16' 
nsplits = 20
dayint = 32   # from control file
ndays = 365   # total number of days across simulation period 

print '   Output path:', path_to_output

## IMPORT GRID-VEG DATA

path_to_gridveg = '/home/563/ajn563/short/fluoro2/control_bethy/'

#f = netCDF4.Dataset('%sbethy_loresgrid3veg_faparratio.nc' % path_to_gridveg,'r')    # lores
f = netCDF4.Dataset('%sbethy_grid3veg_v2.nc' % path_to_gridveg,'r')        # hires

gridnum = f.variables['gridnum'][::-1,:]
vtype = f.variables['type'][:,::-1,:]
frac = f.variables['frac'][:,::-1,:]
lats = f.variables['lat'][::-1]
lons = f.variables['lon'][:]

f.close()


## PULL OUT GRID-VEG DATA

ind = np.where(gridnum != 0)

gridnum_x = gridnum[ind]

frac1_slice = frac[0,:,:]
vtype1_slice = vtype[0,:,:]
frac1 = frac1_slice[ind]
vtype1 = vtype1_slice[ind]

frac2_slice = frac[1,:,:]
vtype2_slice = vtype[1,:,:]
frac2 = frac2_slice[ind]
vtype2 = vtype2_slice[ind]

frac3_slice = frac[2,:,:]
vtype3_slice = vtype[2,:,:]
frac3 = frac3_slice[ind]
vtype3 = vtype3_slice[ind]

frac_gp = np.vstack((frac1,frac2,frac3))
vtype_gp = np.vstack((vtype1,vtype2,vtype3))


## UNRAVEL GRID-VEG DATA INTO VP-length PANDAS ARRAYS

ngp = 3462               # lores = 170; hires = 3462 
nvp = 8776              # lores = 506; hires = 8776
vegp = np.zeros(nvp)
vtype_vp = np.zeros(nvp)
frac_vp = np.zeros(nvp)
gridp_vp = np.zeros(nvp)

jvp = 0   # counter for vegetation points

for i in range(ngp):
    for ifrac in range(3):
        if vtype_gp[ifrac,i] == 0:
            print '   gridp',i+1,'; < 3 veg points; ignoring veg points with 0 fractional coverage'
        else:
            jvp += 1
            vegp[jvp-1] = jvp
            vtype_vp[jvp-1] = vtype_gp[ifrac,i]
            frac_vp[jvp-1] = frac_gp[ifrac,i]
            gridp_vp[jvp-1] = gridnum_x[i]

# Put into pandas dataframe
d = {'vegp': vegp, 'gridp': gridp_vp, 'vtype': vtype_vp, 'frac': frac_vp}
df_vp = pd.DataFrame(d)


## IMPORT BETHY-SCOPE OUTPUT FILES

ntimesteps = 8760

gpp_catblocks = np.zeros((nvp,ntimesteps))    # shape(nvps,ntimesteps across whole sim period)
sif_catblocks = np.zeros((nvp,ntimesteps))

blocks = np.arange(0,nvp,nvp/nsplits)
blocks[-1] = nvp      # Set the last block to end at nvps (i.e. 506 for lores; 8776 for hires)
#print 'Blocks:',blocks

for iblock in range(nsplits):
#    print 'block:',iblock+1
    #Set index bounds for each block
    i1block = blocks[iblock]
    i2block = blocks[iblock+1]
    #Import GPP and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rgppfluo_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    gpp = f.variables['rgppfluo_diurnal'][:]
    f.close()
    vpblock = gpp[i1block:i2block,:]
    gpp_catblocks[i1block:i2block,:] = vpblock
    #Import SIF and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rfluo_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    sif = f.variables['rfluo_diurnal'][:]
    f.close()
    vpblock = sif[i1block:i2block,:]
    sif_catblocks[i1block:i2block,:] = vpblock


## CREATE DATETIME ARRAY

ind = np.arange(0,ndays,dayint)
#ind = np.arange(15,365,30)

base = datetime(2010,1,1,0)
arr = np.array([base + timedelta(days=i) for i in range(ndays)])
#arr = np.array([base + timedelta(days=i) for i in range(365)])

date_arr = []

for i in range(len(ind)):
    date_arr.append(pd.date_range(str(arr[ind][i]),freq='H',periods=24))

date_arr = np.datetime_as_string(date_arr,timezone='UTC')
date_arr = date_arr.flatten().astype('datetime64')
dates_pd = pd.DatetimeIndex(date_arr)

## WRITE TO NETCDF FILE

dataset = netCDF4.Dataset('%spar%s_scope_out_diurnal.nc' % (path_to_output,param), 'w', format='NETCDF4_CLASSIC')

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
nc_times[:] = date2num(dates_pd.astype(datetime), units=nc_times.units, calendar=nc_times.calendar)
nc_vegp[:] = df_vp.vegp.values
nc_gridp[:] = df_vp.gridp.values
nc_frac[:] = df_vp.frac.values
nc_vtype[:] = df_vp.vtype.values
nc_sif[:,:] = sif_catblocks
nc_gpp[:,:] = gpp_catblocks

dataset.close()

