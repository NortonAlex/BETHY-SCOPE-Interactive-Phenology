PRO ipcc_echam5_a1

;dir = '/Users/wolfgang/Models/CCDAS/input/climate/'
dir = '/Volumes/Work/ccdas_input/climate/'
model = 'ECHAM5'
scen = 'A1'
res = 'T63'
out_res = '2x2L'
out_missing = 1e6
sout_missing = '1e6'
name = '_A1.nc'
start_year = 2001
start_year_out = 2008
end_year_out = 2099

;ncdf_read, file = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', vari = 'gridnum', grid
ncdf_read, file = '/Volumes/Work/wolfgang/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', vari = 'gridnum', grid
land = grid.gridnum ne -9999.
oc = grid.gridnum eq -9999.
nx = n_elements(grid.lon)
ny = n_elements(grid.lat)

file_clt = dir + 'scenarios/' + model + '/clt'+name
file_tas = dir + 'scenarios/' + model + '/tas'+name
file_pr = dir + 'scenarios/' + model + '/pr'+name
filen = dir + model + '_' + scen + '_' + $
  strcompress(start_year_out,/rem) + '_' + $
  strcompress(end_year_out,/rem) + '_monthly.nc'

ncdf_read, file = file_clt, attributes=attributes
print,string(attributes.global.source)
print,string(attributes.global.title)
print,attributes.time.units
print,attributes.clt.units

ncdf_read, file = file_clt, data, vari=['clt']
ni = n_elements(data.lon)
nj = n_elements(data.lat)
nm = (end_year_out - start_year_out + 1) * 12
m0 = (start_year_out - start_year) * 12
m1 = (end_year_out - start_year) * 12 + 11
clt = regrid3d(data.clt[*,*,m0:m1], res, out_res)
; convert from cloud percent to cloud fraction
clt = reverse(clt,2) * rebin(land, nx, ny, nm) / 100. $
      + out_missing * rebin(oc, nx, ny, nm)
print, 'done clt'

ncdf_read, file = file_tas, attributes=attributes
print,attributes.tas.units

ncdf_read, file = file_tas, data, vari = ['tas']
tas = regrid3d(data.tas[*,*,m0:m1], res, out_res)
; convert temperature from K to degrees Celsius
tas = reverse (tas - 273.15, 2) * rebin(land, nx, ny, nm)  $
      + out_missing * rebin(oc, nx, ny, nm)
print, 'done tas'

ncdf_read, file = file_pr, attributes=attributes
print,attributes.pr.units

ncdf_read, file = file_pr, data, vari = ['pr']
pr = regrid3d(data.pr[*,*,m0:m1], res, out_res)
; convert precipitation from kg / m^2 / s to mm / day
pr = reverse(pr,2) * rebin(land, nx, ny, nm)  * 3600. * 24. $
     + out_missing * rebin(oc, nx, ny, nm)
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
NCDF_VARPUT, id, vtime, (indgen(nm)+0.5)/12.+start_year_out
NCDF_VARPUT, id, vtemp, tas
NCDF_VARPUT, id, vprecip, pr
NCDF_VARPUT, id, vcloud, clt

   ; close file
NCDF_CLOSE, id

END
