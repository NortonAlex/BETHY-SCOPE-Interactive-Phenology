#! /usr/bin/python

import sys
import numpy as np
import netCDF4
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num

## Command line arguments
#  first argument = experiment name (e.g. par01)
#  second argument = spatial resolution of model run (i.e. lores or hires)
#  third argument = number of blocks/splits 

# Command line arguments 

if len(sys.argv) > 1:
   for arg in sys.argv:
      print 'Command line argument:', arg
else:
   print 'Require parameter number as command line argument... exiting'
   sys.exit() 

# SET SOME VARIABLES

resolution = str(sys.argv[2])    # 'lores' (10x7.5) or 'hires' (2x2)
path_to_output = '/home/563/ajn563/short/fluoro2/output/%s/' % str(sys.argv[1])
exp = str(sys.argv[1])    # 'par16' 
nsplits = int(sys.argv[3])  
base = datetime(2015,1,1,0)    # starting datetime for simulated period
dayint = 1   # if using monthly mean forcing for photosynthesis use dayint=0, otherwise use same as in control file
ndays = 365   # total number of days across simulation period 
path_to_vpblockfile = '/home/563/ajn563/short/fluoro2/control_blockruns/block_vplist.txt'

print '   Output path:', path_to_output

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
    ngp = 170
    nvp = 506
elif resolution == 'hires':
    print '   high-resolution grid scale specified'
    f = netCDF4.Dataset('%sbethy_grid3veg_v2.nc' % path_to_gridveg,'r')        # hires
    gridnum = f.variables['gridnum'][::-1,:]
    vtype = f.variables['type'][:,::-1,:]
    frac = f.variables['frac'][:,::-1,:]
    lats = f.variables['lat'][::-1]
    lons = f.variables['lon'][:]
    f.close()
    ngp = 3462
    nvp = 8776 
else:
    print '   No or incorrect model grid resolution specified. Must specify as 2nd command line argument on call... exiting'
    sys.exit()

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

ntimesteps = ndays*24    # number of hourly time-steps across simulation period

# Initiate output arrays
gpp_catblocks = np.zeros((nvp,ntimesteps))    # shape(nvps,ntimesteps across whole sim period)
sif_catblocks = np.zeros((nvp,ntimesteps))
lai_catblocks = np.zeros((nvp,ntimesteps))
apar_catblocks = np.zeros((nvp,ntimesteps))
aparcab_catblocks = np.zeros((nvp,ntimesteps))
par_catblocks = np.zeros((nvp,ntimesteps))
swdown_catblocks = np.zeros((nvp,ntimesteps))
pres_catblocks = np.zeros((nvp,ntimesteps))
ta_catblocks = np.zeros((nvp,ntimesteps))
ea_catblocks = np.zeros((nvp,ntimesteps))

# Fill initiated output arrays with missing-value
gpp_catblocks.fill(-999.9)
sif_catblocks.fill(-999.9)
lai_catblocks.fill(-999.9)
apar_catblocks.fill(-999.9)
aparcab_catblocks.fill(-999.9)
par_catblocks.fill(-999.9)
swdown_catblocks.fill(-999.9)
pres_catblocks.fill(-999.9)
ta_catblocks.fill(-999.9)
ea_catblocks.fill(-999.9)

#blocks = np.arange(0,nvp,nvp/nsplits)
#blocks[-1] = nvp      # Set the last block to end at nvps (i.e. 506 for lores; 8776 for hires)
#print 'Blocks:',blocks

# Read in selected veg-points 
f = open('%s' % path_to_vpblockfile)
lines = f.read().splitlines()
iline = 1

ivp_all = np.zeros(0)

for iblock in range(nsplits):
#    print 'block:',iblock+1
    #Set index bounds for each block
    xline = lines[iline].split(' ')[0:-1]    # exclude index point -1 as it is a space ' '
    vp = np.array([int(i) for i in xline],dtype=int)    # selected veg-points that were actually run in scope
    print '   veg-points:',vp
    ivp = vp - 1    # index points of selected veg-points
    ivp_all = np.append(ivp_all,ivp)
    #Import GPP and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rgppfluo_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    gpp = f.variables['rgppfluo_diurnal'][:]
    f.close()
    vpblock = gpp[ivp,:]
    gpp_catblocks[ivp,:] = vpblock
    #Import SIF and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rfluo_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    sif = f.variables['rfluo_diurnal'][:]
    f.close()
    vpblock = sif[ivp,:]
    sif_catblocks[ivp,:] = vpblock
    #Import LAI and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rlai_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    lai = f.variables['rlai_diurnal'][:]
    f.close()
    vpblock = lai[ivp,:]
    lai_catblocks[ivp,:] = vpblock
    #Import SCOPE APAR and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rapar_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    apar = f.variables['rapar_diurnal'][:]
    f.close()
    vpblock = apar[ivp,:]
    apar_catblocks[ivp,:] = vpblock
    #Import SCOPE APAR (Cab) and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/raparcab_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    aparcab = f.variables['raparcab_diurnal'][:]
    f.close()
    vpblock = aparcab[ivp,:]
    aparcab_catblocks[ivp,:] = vpblock
    #Import SCOPE PAR (incident par) and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rpar_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    par = f.variables['rpar_diurnal'][:]
    f.close()
    vpblock = par[ivp,:]
    par_catblocks[ivp,:] = vpblock
    #Import SCOPE SWDOWN (shortwave radiation down) and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rswdown_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    swdown = f.variables['rswdown_diurnal'][:]
    f.close()
    vpblock = swdown[ivp,:]
    swdown_catblocks[ivp,:] = vpblock
    #Import SCOPE PRES (PA, air pressure, hPa) and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rpres_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    pres = f.variables['rpres_diurnal'][:]
    f.close()
    vpblock = pres[ivp,:]
    pres_catblocks[ivp,:] = vpblock
    #Import SCOPE TA (air temperature, oC) and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rta_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    ta = f.variables['rta_diurnal'][:]
    f.close()
    vpblock = ta[ivp,:]
    ta_catblocks[ivp,:] = vpblock
    #Import SCOPE EA (vapour pressure, hPa) and splice into concatenated array
    f = netCDF4.Dataset('%sjob{0:02d}/rea_diurnal.nc'.format(iblock+1) % path_to_output,'r')
    ea = f.variables['rea_diurnal'][:]
    f.close()
    vpblock = ea[ivp,:]
    ea_catblocks[ivp,:] = vpblock
    #
    iline+=2

# Only keep the time slices where data exists. Missing value = -999 (see code).
ivp_first = int(ivp_all[0])    # first veg-point with data (i.e. run in selected vp simulation)
inds = np.squeeze(np.where(gpp_catblocks[ivp_first,:] > -990))    # indexes of time slices with real values (i.e. >-990) to keep.
#inds = np.array([13,757,1429,2173,2893,3637,4357,5101,5845,6565,7309,8029])

if np.size(inds) == 0: 
    print 'Error: no valid scope values at any time-step' 
print 'tinds: size,min,max ',np.size(inds),np.min(inds),np.max(inds)
print '  sif_catblocks: min, max:',np.nanmin(sif_catblocks),np.nanmax(sif_catblocks)
print '     size of non-zero values:',np.shape(np.where(sif_catblocks > -990))

gpp_out = gpp_catblocks[:,inds]
sif_out = sif_catblocks[:,inds]
lai_out = lai_catblocks[:,inds]
apar_out = apar_catblocks[:,inds]
aparcab_out = aparcab_catblocks[:,inds]
par_out = par_catblocks[:,inds]
swdown_out = swdown_catblocks[:,inds]
pres_out = pres_catblocks[:,inds]
ta_out = ta_catblocks[:,inds]
ea_out = ea_catblocks[:,inds]

## CREATE DATETIME ARRAY

if dayint == 0:
    print 'dayint = 0'
    # Assign the DOY index for the 1st of every month
    first_idom = np.array([0,31,59,90,120,151,181,212,243,273,304,334])

    print '    hello, I hope these runs used monthly mean forcing for GPP'
    print '    ...setting the datetime array to be first day of each month'
    arr = np.array([base + timedelta(days=i) for i in range(ndays)])

    date_arr = []

    for i in range(len(first_idom)):
        date_arr.append(pd.date_range(str(arr[first_idom][i]),freq='H',periods=24))

else:
    ind = np.arange(0,ndays,dayint)
    arr = np.array([base + timedelta(days=i) for i in range(ndays)])
    date_arr = []

    for i in range(len(ind)):
        date_arr.append(pd.date_range(str(arr[ind][i]),freq='H',periods=24))

date_arr = np.datetime_as_string(date_arr,timezone='UTC')
date_arr = date_arr.flatten().astype('datetime64')
dates_pd = pd.DatetimeIndex(date_arr)


print 'dates_pd: size,0,-1 ',np.size(dates_pd),dates_pd[0],dates_pd[-1]
## WRITE TO NETCDF FILE

dataset = netCDF4.Dataset('%s%s_scope_out.nc' % (path_to_output,exp), 'w', format='NETCDF4_CLASSIC')

# Create Dimensions
dataset.createDimension('time',None)
dataset.createDimension('vegp',nvp)

# Create Variables
nc_times = dataset.createVariable('time', np.float32, ('time',)) 
nc_vegp = dataset.createVariable('vegp', np.int32, ('vegp',),fill_value=-999.9)
nc_gridp = dataset.createVariable('gridp', np.int32, ('vegp',),fill_value=-999.9)
nc_frac = dataset.createVariable('frac', np.float32, ('vegp',),fill_value=-999.9)
nc_vtype = dataset.createVariable('vtype', np.int32, ('vegp',),fill_value=-999.9)
nc_sif = dataset.createVariable('SIF', np.float32, ('time','vegp'),fill_value=-999.9)
nc_gpp = dataset.createVariable('GPP', np.float32, ('time','vegp'),fill_value=-999.9)
nc_lai = dataset.createVariable('LAI', np.float32, ('time','vegp'),fill_value=-999.9)
nc_apar = dataset.createVariable('APAR', np.float32, ('time','vegp'),fill_value=-999.9)
nc_aparcab = dataset.createVariable('APAR_Cab', np.float32, ('time','vegp'),fill_value=-999.9)
nc_par = dataset.createVariable('PAR', np.float32, ('time','vegp'),fill_value=-999.9)
nc_swdown = dataset.createVariable('SWDOWN', np.float32, ('time','vegp'),fill_value=-999.9)
nc_pres = dataset.createVariable('PRES', np.float32, ('time','vegp'),fill_value=-999.9)
nc_ta = dataset.createVariable('TA', np.float32, ('time','vegp'),fill_value=-999.9)
nc_ea = dataset.createVariable('EA', np.float32, ('time','vegp'),fill_value=-999.9)

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
nc_lai.long_name = 'Leaf-Area Index'
nc_lai.units = 'm2 m-2'
nc_apar.long_name = 'Canopy Absorbed Photosynthetically Active Radiation'
nc_apar.units = 'umol photons m-2 s-1'
nc_aparcab.long_name = 'Canopy Absorbed Photosynthetically Active Radiation by Chlorophyl a-b'
nc_aparcab.units = 'umol photons m-2 s-1' 
nc_par.long_name = 'Incident Photosynthetically Active Radiation above Canopy'
nc_par.units = 'umol photons m-2 s-1'
nc_swdown.long_name = 'Shortwave Down Radiation'
nc_swdown.units = 'W m-2'
nc_pres.long_name = 'Air Pressure'
nc_pres.units = 'hPa'
nc_ta.long_name = 'Air Temperature'
nc_ta.units = 'oC'
nc_ea.long_name = 'Vapour Pressure'
nc_ea.units = 'hPa'

# Fill variables with data
nc_times[:] = date2num(dates_pd[inds].astype(datetime), units=nc_times.units, calendar=nc_times.calendar)
nc_vegp[:] = df_vp.vegp.values
nc_gridp[:] = df_vp.gridp.values
nc_frac[:] = df_vp.frac.values
nc_vtype[:] = df_vp.vtype.values
nc_sif[:,:] = np.swapaxes(sif_out,0,1)
nc_gpp[:,:] = np.swapaxes(gpp_out,0,1)
nc_lai[:,:] = np.swapaxes(lai_out,0,1)
nc_apar[:,:] = np.swapaxes(apar_out,0,1)
nc_aparcab[:,:] = np.swapaxes(aparcab_out,0,1)
nc_par[:,:] = np.swapaxes(par_out,0,1)
nc_swdown[:,:] = np.swapaxes(swdown_out,0,1)
nc_pres[:,:] = np.swapaxes(pres_out,0,1)
nc_ta[:,:] = np.swapaxes(ta_out,0,1)
nc_ea[:,:] = np.swapaxes(ea_out,0,1)

dataset.close()

