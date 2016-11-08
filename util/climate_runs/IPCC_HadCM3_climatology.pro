PRO ipcc_hadcm3_climatology

;device,decompose=0
;loadct,41

;dir = '/Users/wolfgang/Models/CCDAS/input/climate/'
dir = '/Volumes/Work/ccdas_input/climate/'
model = 'HadCM3'
scen = 'A1'
res = '3.75x1.25'
out_res = '2x2L'
out_missing = 1e6
sout_missing = '1e6'
name1 = '_A1_20th_century.nc'
start_year1 = 1860
name2 = '_A1.nc'
start_year2 = 2000
year0 = 1979 ; first year of climatology
year1 = 2007 ; last year of climatology

;ncdf_read, file = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', vari = 'gridnum', grid
ncdf_read, file = '/Volumes/Work/wolfgang/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', vari = 'gridnum', grid
land = grid.gridnum ne -9999.
oc = grid.gridnum eq -9999.
nx = n_elements(grid.lon)
ny = n_elements(grid.lat)

file_clt1 = dir + 'scenarios/' + model + '/clt'+name1
file_clt2 = dir + 'scenarios/' + model + '/clt'+name2
file_tas1 = dir + 'scenarios/' + model + '/tas'+name1
file_tas2 = dir + 'scenarios/' + model + '/tas'+name2
file_pr1 = dir + 'scenarios/' + model + '/pr'+name1
file_pr2 = dir + 'scenarios/' + model + '/pr'+name2
filen = dir + model + '_' + scen + '_'+strcompress(year0,/rem) $
         +'_'+strcompress(year1,/rem)+'_monthly_climatology.nc'

ncdf_read, file = file_clt1, attributes=attributes
print,string(attributes.global.source)
print,attributes.time.units
print,attributes.clt.units

ncdf_read, file = file_clt2, attributes=attributes
print,string(attributes.global.source)
print,string(attributes.global.comment)
print,attributes.time.units
print,attributes.clt.units

ncdf_read, file = file_clt1, data1, vari=['clt']
ncdf_read, file = file_clt2, data2, vari=['clt']
ni = n_elements(data1.lon)
nj = n_elements(data1.lat)
nm1 = n_elements(data1.time)
nm2 = n_elements(data2.time)
nm = (year1 - year0 + 1) * 12
m1 = (year0-start_year1)*12     ; first month used for climatology
m2 = (year1-start_year2)*12 + 11 ; last month used for climatology (in other data set)
clt0 = fltarr(ni,nj,nm)
clt0[*,*,0:nm1-m1-1] = data1.clt[*,*,m1:nm1-1]
clt0[*,*,nm1-m1:nm-1] = data2.clt[*,*,0:m2]
clt0 = reform(clt0,ni,nj,12,year1-year0+1)
clt0 = total (clt0, 4) / (year1-year0+1)
clt = fltarr (ni, 2*nj-2, 12)
clt[*,0,*] = clt0[*,0,*]
clt[*,1:2*nj-4,*] = rebin(clt0[*,1:nj-2,*],ni,2*nj-4,12)
clt[*,2*nj-3,*] = clt0[*,nj-1,*]
clt = regrid3d(clt, res, out_res)
; convert from cloud percent to cloud fraction
clt = clt * rebin(land, nx, ny, 12) / 100. $
      + out_missing * rebin(oc, nx, ny, 12)
print, 'done clt'

ncdf_read, file = file_tas1, attributes=attributes
print,attributes.tas.units

ncdf_read, file = file_tas1, data1, vari = ['tas']
ncdf_read, file = file_tas2, data2, vari = ['tas']
tas0 = fltarr(ni,nj,nm)
tas0[*,*,0:nm1-m1-1] = data1.tas[*,*,m1:nm1-1]
tas0[*,*,nm1-m1:nm-1] = data2.tas[*,*,0:m2]
tas0 = reform(tas0,ni,nj,12,year1-year0+1)
; check: 
;plot,total(total(total(tas0,1),1),1)/12/ni/nj-273.15,/ynozero
tas0 = total (tas0, 4) / (year1-year0+1)
tas = fltarr (ni, 2*nj-2, 12)
tas[*,0,*] = tas0[*,0,*]
tas[*,1:2*nj-4,*] = rebin(tas0[*,1:nj-2,*],ni,2*nj-4,12)
tas[*,2*nj-3,*] = tas0[*,nj-1,*]
tas = regrid3d(tas, res, out_res)
; checking:
;tvf,/bar,min=-20,max=30,fo='(I4)',n=6,-273.15+tas[*,*,6]
;tvf,/bar,min=-20,max=30,fo='(I4)',n=6,-273.15+tas0[*,*,6]
; convert temperature from K to degrees Celsius
tas = (tas - 273.15) * rebin(land, nx, ny, 12)  $
      + out_missing * rebin(oc, nx, ny, 12)
; checking
;tvf,/map,mis=1e6,min=-20,max=30,/bar,n=6,fo='(I4)',tas[*,*,0]
print, 'done tas'

ncdf_read, file = file_pr1, attributes=attributes
print,attributes.pr.units

ncdf_read, file = file_pr1, data1, vari = ['pr']
ncdf_read, file = file_pr2, data2, vari = ['pr']
pr0 = fltarr(ni,nj,nm)
pr0[*,*,0:nm1-m1-1] = data1.pr[*,*,m1:nm1-1]
pr0[*,*,nm1-m1:nm-1] = data2.pr[*,*,0:m2]
pr0 = reform(pr0,ni,nj,12,year1-year0+1)
pr0 = total (pr0, 4) / (year1-year0+1)
pr = fltarr (ni, 2*nj-2, 12)
pr[*,0,*] = pr0[*,0,*]
pr[*,1:2*nj-4,*] = rebin(pr0[*,1:nj-2,*],ni,2*nj-4,12)
pr[*,2*nj-3,*] = pr0[*,nj-1,*]
pr = regrid3d(pr, res, out_res)
; convert precipitation from kg / m^2 / s to mm / day
pr = pr * rebin(land, nx, ny, 12)  * 3600. * 24. $
     + out_missing * rebin(oc, nx, ny, 12)
print, 'done pr'

id = NCDF_CREATE(filen, /CLOBBER)

  ; define dimensions
dlon   = NCDF_DIMDEF(id, 'lon', nx)
dlat   = NCDF_DIMDEF(id, 'lat', ny)
dtime  = NCDF_DIMDEF(id, 'time', 12)

vlon = NCDF_VARDEF(id, 'lon', [dlon], /FLOAT)
NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

vlat = NCDF_VARDEF(id, 'lat', [dlat], /FLOAT)
NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
NCDF_ATTPUT, id, vtime, 'long_name', 'month_of_year', /CHAR
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
