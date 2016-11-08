PRO VIC2lores

yr0 = 1979
yr1 = 2008
s_yrs = strcompress(yr0,/rem)+'-'+strcompress(yr1,/rem)
path = '/Volumes/Work/data/VIC_drivers/revision_nov_2008/'
names = ['cloud','precip','swdown','tmax','tmin','lwdown','vapor']
units = ['fraction','mm/day', 'W/m2', 'deg C', 'deg C','W/m2', 'kPa']
long_names = ['cloud cover', 'precipitation total', 'incoming shortwave radiation', 'daily maximum 2m temperature', 'daily minimum 2m temperature','incoming longwave radiation', 'daily atmospheric saturated vapor pressure']
missing = 1e20

; get the land-sea mask and lon/lat coordinates in TM2 resolution (="lo res")
ncdf_read, file='/Users/wolfgang/Models/CCDAS/ppm/control_bethy/adbethy_loresgrid3veg_79-00.nc', variables=['gridnum','lon','lat'], lores_in
lsm_lores = lores_in.gridnum ne -9999.

; loop over different variables
for cnt=0,n_elements(names)-1 do begin
;for cnt=0,0 do begin

; read various fields from the 2deg by 2deg input file in NetCDF
fname = path+'forcing.2x2.'+names[cnt]+'.'+s_yrs+'.nc'
print, 'Reading file ', fname, '...'
ncdf_read, file=fname, variables = names[cnt], /no_struct, in
ncdf_read, file=fname, variables = 'lon', /no_struct, lon
ncdf_read, file=fname, variables = 'lat', /no_struct, lat
ncdf_read, file=fname, variables = 'time', /no_struct, time

; dimensions of input and output fields
s = size (in)
s_lores = size (lsm_lores)
; mask = rebin (lsm_lores, s_lores[1], s_lores[2], s[3])

; convert to lo res
; out = bethy2tm (in, miss=missing) ; alternative formulation - faster, but needs more memory
; out = out * mask + missing * (mask eq 0)
out = fltarr (s_lores[1], s_lores[2], s[3])
for t = 0, s[3]-1 do out[*,*,t] = bethy2tm (in[*,*,t], miss=missing) * lsm_lores + missing * (lsm_lores eq 0)

; display sample of result
; device,decompose=0
; loadct, 13
; window, 0
; tvf,in[*,*,s[3]-1],/bar,missing=missing
; window, 1
; tvf,out[*,*,s[3]-1],/bar,missing=missing

; create the NetCDF output file with new resolution
fnew = path+'forcing.lores.'+names[cnt]+'.'+s_yrs+'.nc'
id = NCDF_CREATE(fnew, /CLOBBER)

; define dimensions
dlon   = NCDF_DIMDEF(id, 'lon', s_lores[1])
dlat   = NCDF_DIMDEF(id, 'lat', s_lores[2])
dtime  = NCDF_DIMDEF(id, 'time', s[3])

; description of first dimension
vlon = NCDF_VARDEF(id, 'lon', [dlon], /FLOAT)
NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

; description of second dimension
vlat = NCDF_VARDEF(id, 'lat', [dlat], /FLOAT)
NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

; description of third dimension
vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
NCDF_ATTPUT, id, vtime, 'long_name', 'time', /CHAR
NCDF_ATTPUT, id, vtime, 'units', 'days since 1-1-1979', /CHAR

; define variable
vvar = NCDF_VARDEF(id, names[cnt], [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vvar, 'long_name', long_names[cnt], /CHAR
NCDF_ATTPUT, id, vvar, 'units', units[cnt], /CHAR
NCDF_ATTPUT, id, vvar, 'missing_value', '1e20', /FLOAT

; change to DATA mode
NCDF_CONTROL, id, /ENDEF

; put data into file
NCDF_VARPUT, id, vlon, lores_in.lon
NCDF_VARPUT, id, vlat, lores_in.lat
NCDF_VARPUT, id, vtime, time
NCDF_VARPUT, id, vvar, out

; close file
NCDF_CLOSE, id
print, 'File written: ', fnew
print, ' '

; end of loop over files
ENDFOR

; testing the result:
; cnt = 0
; ncdf_read, file=path+'forcing.lores.'+names[cnt]+'.'+s_yrs+'.nc', /all, in_test, $
; attributes=attributes
; tvf,in_test.cloud[*,*,0],/bar,missing=missing
; tvf,in_test.cloud[*,*,s[3]-1],/bar,missing=missing
; help,/str,attributes.lon
; help,/str,attributes.lat
; help,/str,attributes.cloud

END
