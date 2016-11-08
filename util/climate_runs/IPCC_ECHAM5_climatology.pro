PRO ipcc_echam5_climatology

;dir = '/Users/wolfgang/Models/CCDAS/input/climate/'
dir = '/Volumes/Work/ccdas_input/climate/'
model = 'ECHAM5'
scen = 'A1'
res = 'T63'
out_res = '2x2L'
out_missing = 1e6
sout_missing = '1e6'
start_year = 1860
year0 = 1979 ; first year of climatology
year1 = 2007 ; last year of climatology

;ncdf_read, file = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', vari = 'gridnum', grid
ncdf_read, file = '/Volumes/Work/wolfgang/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', vari = 'gridnum', grid
land = grid.gridnum ne -9999.
oc = grid.gridnum eq -9999.
nx = n_elements(grid.lon)
ny = n_elements(grid.lat)

file_clt = dir + 'scenarios/' + model + '/clt_'+scen+'_20th_century.nc'
file_tas = dir + 'scenarios/' + model + '/tas_'+scen+'_20th_century.nc'
file_pr = dir + 'scenarios/' + model + '/pr_'+scen+'_20th_century.nc'
filen = dir + model + '_' + scen + '_'+strcompress(year0,/rem) $
         +'_'+strcompress(year1,/rem)+'_monthly_climatology.nc'

ncdf_read, file = file_clt, attributes=attributes
print,string(attributes.global.source)
print,string(attributes.global.title)
print,attributes.time.units
print,attributes.clt.units

ncdf_read, file = file_clt, data, vari=['clt']
ni = n_elements(data.lon)
nj = n_elements(data.lat)

m0 = (year0-start_year)*12     ; first month used for climatology
m1 = (year1-start_year)*12 + 11 ; last month used for climatology
clt = reform(data.clt[*,*,m0:m1],ni,nj,12,year1-year0+1)
clt = total (clt, 4) / (year1-year0+1)
clt = regrid3d(clt, res, out_res)
; convert from cloud percent to cloud fraction
clt = reverse(clt,2) * rebin(land, nx, ny, 12) / 100. $
      + out_missing * rebin(oc, nx, ny, 12)
print, 'done clt'

ncdf_read, file = file_tas, attributes=attributes
print,attributes.tas.units

ncdf_read, file = file_tas, data, vari = ['tas']
tas = reform(data.tas[*,*,m0:m1],ni,nj,12,year1-year0+1)
tas = total (tas, 4) / (year1-year0+1)
tas = regrid3d(tas, res, out_res)
; convert temperature from K to degrees Celsius
tas = reverse (tas - 273.15, 2) * rebin(land, nx, ny, 12)  $
      + out_missing * rebin(oc, nx, ny, 12)
print, 'done tas'

ncdf_read, file = file_pr, attributes=attributes
print,attributes.pr.units

ncdf_read, file = file_pr, data, vari = ['pr']
pr = reform(data.pr[*,*,m0:m1],ni,nj,12,year1-year0+1)
pr = total (pr, 4) / (year1-year0+1)
pr = regrid3d(pr, res, out_res)
; convert precipitation from kg / m^2 / s to mm / day
pr = reverse(pr,2) * rebin(land, nx, ny, 12)  * 3600. * 24. $
     + out_missing * rebin(oc, nx, ny, 12)
print, 'done pr'

id = NCDF_CREATE(filen, /CLOBBER)

  ; define dimensions
dlon   = NCDF_DIMDEF(id, 'lon', nx)
dlat   = NCDF_DIMDEF(id, 'lat', ny)
dtime  = NCDF_DIMDEF(id, 'time', /UNLIMITED)

vlon = NCDF_VARDEF(id, 'lon', [dlon], /FLOAT)
NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

vlat = NCDF_VARDEF(id, 'lat', [dlat], /FLOAT)
NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
NCDF_ATTPUT, id, vtime, 'long_name', 'time since Jan 2000', /CHAR
NCDF_ATTPUT, id, vtime, 'units', 'month', /CHAR

   ; define variables
  
vtemp = NCDF_VARDEF(id, 'temp', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vtemp, 'long_name', 'Mean monthly surface air temperature', /CHAR
NCDF_ATTPUT, id, vtemp, 'units', 'degrees C', /CHAR
NCDF_ATTPUT, id, vtemp, 'missing_value', sout_missing, /FLOAT

vprecip = NCDF_VARDEF(id, 'precip', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vprecip, 'long_name', 'Mean monthly precipitation', /CHAR
NCDF_ATTPUT, id, vprecip, 'units', 'mm/day', /CHAR
NCDF_ATTPUT, id, vprecip, 'missing_value', sout_missing, /FLOAT

vcloud = NCDF_VARDEF(id, 'cloud', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vcloud, 'long_name', 'Mean monthly total cloud fraction', /CHAR
NCDF_ATTPUT, id, vcloud, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vcloud, 'missing_value', sout_missing, /FLOAT

   ; define global attributes
NCDF_ATTPUT, id, /GLOBAL, 'title', model + ' ' + scen + ' regridded on 2x2 deg', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'author', 'Marko Scholze & Wolfgang Knorr', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR

   ; change to DATA mode and put data into file
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, vlon, grid.lon
NCDF_VARPUT, id, vlat, grid.lat
NCDF_VARPUT, id, vtime, indgen(12) + 1
NCDF_VARPUT, id, vtemp, tas
NCDF_VARPUT, id, vprecip, pr
NCDF_VARPUT, id, vcloud, clt

   ; close file
NCDF_CLOSE, id

END
