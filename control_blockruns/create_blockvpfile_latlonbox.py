#! /usr/bin/python

import numpy as np
import netCDF4
import pandas as pd
import sys

## Read me:
# This script produces a text file that is used to run only selected region of interest
# in the BETHY-SCOPE. This text file is then read into the code and fluo is run on these 
# veg-points only. In the text file:
# - Every odd line is the number of veg-points per block
# - Every even line is the veg-points for that block

print 'Writing vp-block file with selected lat-lon box region'

## Command line arguments:
## first argument = resolution (i.e. lores or hires)
## second argument = nsplits  (i.e. nblocks)


if len(sys.argv) > 1:
   for arg in sys.argv:
      print 'Command line argument:', arg
else:
   print 'Require multiple command line arguments... exiting'
   sys.exit()

resolution = str(sys.argv[1])
nsplits = int(sys.argv[2])

#pft_list = sys.argv[3].split(',')    # third argument
#pft_sel = np.array([int(i) for i in pft_list])    # convert string list into int array
## Select region by a lat-lon box
lat0 = 46
lat1 = 50
lon0 = -100
lon1 = -66

print '   Model resolution selected:',resolution
#print '             PFT(s) selected:',pft_sel
print '                     nsplits:',nsplits
print '   Selected region,  Lat-Lon Box:'
print '              lats',lat0,'to',lat1
print '              lons',lon0,'to',lon1

## Grid file info

if resolution == 'lores':
    print 'Low-resolution (10x7.5) model scale specified'
    gridfile = '/home/563/ajn563/short/fluoro2/control_bethy/bethy_loresgrid3veg_faparratio.nc'
    nvp = 506
    f = netCDF4.Dataset('%s' % gridfile,'r')
    gridnum = f.variables['gridnum'][::-1,:]
    vtype = f.variables['type'][:,::-1,:]
    frac = f.variables['frac'][:,::-1,:]
    vnum = f.variables['vnum'][:]
    f.close()
elif resolution == 'hires':
    print 'High-resolution (2x2) model scale specified'
    gridfile = '/home/563/ajn563/short/fluoro2/control_bethy/bethy_grid3veg_v2.nc'
    nvp = 8776
    f = netCDF4.Dataset('%s' % gridfile,'r')
    gridnum = f.variables['gridnum'][::-1,:]
    vtype = f.variables['type'][:,::-1,:]
    frac = f.variables['frac'][:,::-1,:]
    vnum = f.variables['vnum'][:]
    f.close()
else:
    print 'No resolution specified... exiting'
    sys.exit()

## - Get land-ocean grid mask

frac_nan = np.ma.getmask(np.isnan(frac[0,:,:]))
print 'shape:', np.shape(frac_nan)
tlen = 12
frac_stacked = np.array([frac_nan for _ in range(tlen)])
print 'shape:', np.shape(frac_stacked)
oceanmask = frac_stacked
igp_landpoints = np.where(oceanmask == False)
igp_landpoints_geog = np.where(oceanmask[0,:,:] == False)

## - Get lats and lons per grid-point
path = '/home/563/ajn563/short/fluoro2/control_blockruns/'

f = netCDF4.Dataset('%spar4000_rfluo_stitched.nc' % path,'r')

lats = f.variables['latitude'][::-1]
lons = f.variables['longitude'][:]

#print 'Before roll:'
#print 'lats:', np.shape(lats)
#print 'lons:', np.shape(lons)
#print '   lats from',lats[0],' to ',lats[-1]
#print '   lons from',lons[0],' to ',lons[-1]

f.close()

latstep = np.abs(lats[1] - lats[0])
lonstep = np.abs(lons[1] - lons[0])

maplats = lats
maplons = lons

xlon, ylat = np.meshgrid(lons, lats)

#igp_landpoints = np.where(oceanmask == False)
#igp_landpoints_geog = np.where(oceanmask[0,:,:] == False)
#print 'shape of mask:', np.shape(igp_landpoints)

lats_gps = ylat[igp_landpoints_geog]
lons_gps = xlon[igp_landpoints_geog]

## - Get 

ind = np.where(gridnum != 0)    # land mask (selects only land grid points)

gridnum = gridnum[ind]

print '  shape of gridnum:',np.shape(gridnum),' first and last elements:',gridnum[0],gridnum[-1]

ngp = len(gridnum)

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
            gridp_vp[jvp-1] = gridnum[i]

# Put into pandas dataframe
d = {'vegp': vegp, 'gridp': gridp_vp, 'vtype': vtype_vp, 'frac': frac_vp}
df_vp = pd.DataFrame(d)

## Select PFT(s) of interest
#if (len(pft_sel) == 1) & (pft_sel[0] == -1):
#    # Select all veg-points
#    print 'Using all veg-points'
#    inds = np.where(df_vp['vtype'].values != 0)
#else:
#    print 'Using PFT(s):',pft_sel
#    inds = np.in1d(df_vp['vtype'].values, pft_sel)

region_inds = np.squeeze(np.where((lats_gps > lat0) & (lats_gps < lat1) & (lons_gps > lon0) & (lons_gps < lon1)))

gridp_sel = gridnum[region_inds]

#gridp_sel = np.unique(df_vp['gridp'].values[region_inds])

print 'first grid point:',gridp_sel[0],' last:',gridp_sel[-1]

# Split up selected grid-points to veg-points then per block
# number of grid-points selected
ngp_sel = np.size(gridp_sel)
print 'number of grid points in the selected region:',ngp_sel

# Set index values for blocks over selected grid-points
#  - If the number of splits it more than the number of selected grid-points (ngp_sel), reset the number of splits to ngp_sel
if ngp_sel < nsplits:
    # one gridp per block
    print 'nsplits is > n grid-points to run'
    nsplits = ngp_sel-1
    iblocksplit = np.arange(0,ngp_sel)
    iblocksplit[-1] = ngp_sel
else:
    iblocksplit = np.array(np.around(np.arange(0,ngp_sel,ngp_sel/float(nsplits))),dtype=int)
    iblocksplit = np.append(iblocksplit,ngp_sel)    # set last block to finish at the last selected point, not at an interval set  by arange

print ''

## Loop over each block, select veg-points and write to text file
outfile = open('block_vplist.txt','w')

for i in range(nsplits):
    print 'Block:',i+1
    igp0 = iblocksplit[i]
    igp1 = iblocksplit[i+1]
    #print '   First gridp in this block:',gridp_sel[igp0]
    #print '   Last  gridp in this block:',gridp_sel[igp1-1]
    #print '   n gridps:',len(gridp_sel[igp0:igp1])
    #print gridp_sel[igp0:igp1]
    gridp_block = gridp_sel[igp0:igp1]
    # Selected veg-points in this block corresponding to the above grid-points
    inds = np.in1d(df_vp['gridp'].values, gridp_block)
    vp_block = df_vp['vegp'].values[inds]
    #print '     Veg points:',vp_block
    ## Save vp_block info in a text file
    ## Odd rows are the number of vp's per block
    outfile.write('%i\n' % len(vp_block))
    ## Even rows are the selected vp's per block
    for j in range(len(vp_block)):
        outfile.write('%i ' % vp_block[j])
    outfile.write('\n')
    
outfile.close()


