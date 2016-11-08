PRO VIC2hires

; take daily forcing data from VIC model and mask out the CCDAS
; hires mask
yr0 = 1979
yr1 = 2008
s_yrs = strcompress(yr0,/rem)+'-'+strcompress(yr1,/rem)
path = '/Volumes/Work/data/VIC_drivers/revision_nov_2008/'
names = ['cloud','precip','swdown','tmax','tmin','lwdown','vapor']
units = ['fraction','mm/day', 'W/m2', 'deg C', 'deg C','W/m2', 'kPa']
long_names = ['cloud cover', 'precipitation total', 'incoming shortwave radiation', 'daily maximum 2m temperature', 'daily minimum 2m temperature','incoming longwave radiation', 'daily atmospheric saturated vapor pressure']
missing = 1e20

; get the land-sea mask and lon/lat coordinates
ncdf_read, file='/Users/wolfgang/Models/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', variables=['gridnum','lon','lat'], hires_in
lsm_hires = hires_in.gridnum ne -9999.

print, 'hires land mask has ', FIX(TOTAL(lsm_hires)), ' grid points.'

; loop over different variables
for cnt=0,n_elements(names)-1 do begin
;for cnt=6,6 do begin

fname = path+'forcing.2x2.'+names[cnt]+'.'+s_yrs+'.nc'
print, 'Reading file ', fname, '...'
; read various fields from the 2deg by 2deg input file in NetCDF
ncdf_read, file=fname, variables = names[cnt], /no_struct, in
ncdf_read, file=fname, variables = 'lon', /no_struct, lon
ncdf_read, file=fname, variables = 'lat', /no_struct, lat
ncdf_read, file=fname, variables = 'time', /no_struct, time

; land mask of input data
lsm_VIC = in[*,*,0] NE missing
IF cnt GT 0 THEN BEGIN
  ; check if mask did not change
  lsm_VIC_new = in[*,*,0] NE missing
  IF total(lsm_VIC_new NE lsm_VIC) GT 0 THEN BEGIN
    print, 'VIC land masks don not match. Stop.'
    stop
  ENDIF
ENDIF ELSE BEGIN
  print, 'VIC land mask has ', FIX(TOTAL(lsm_VIC)), ' grid points.'
ENDELSE

; analyse mismatch
missing_hires = lsm_VIC   GT lsm_hires
missing_VIC   = lsm_hires GT lsm_VIC
print, 'Number of points in VIC but not in hires: ', FIX (TOTAL (missing_hires))
print, 'Number of points in hires but not in VIC: ', FIX (TOTAL (missing_VIC))
IF total(missing_VIC) GT 0 THEN BEGIN
  print, 'Data missing on hires land points. Stop.'
  stop
ENDIF

; applying mask by deleting points from VIC data
out = in
s = size (in)
missing_hires = rebin (missing_hires,s[1],s[2],s[3])
l = where (missing_hires)
out[l] = missing

; testing
lsm_new = out[*,*,0] ne missing
print, 'Points not matching after masking out:      ', fix(total(lsm_new NE lsm_hires))

; construct full time information

; create the NetCDF output file with new land mask

fnew = path+'forcing.hires.'+names[cnt]+'.'+s_yrs+'.nc'
id = NCDF_CREATE(fnew, /CLOBBER)

; define dimensions
dlon   = NCDF_DIMDEF(id, 'lon', s[1])
dlat   = NCDF_DIMDEF(id, 'lat', s[2])
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
NCDF_VARPUT, id, vlon, hires_in.lon
NCDF_VARPUT, id, vlat, hires_in.lat
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
; ncdf_read, file=path+'forcing.hires.'+names[cnt]+'.'+s_yrs+'.nc', /all, in_test, $
; attributes=attributes
; tvf,in_test.cloud[*,*,0],/bar,missing=missing
; tvf,in_test.cloud[*,*,s[3]-1],/bar,missing=missing
; help,/str,attributes.lon
; help,/str,attributes.lat
; help,/str,attributes.cloud

END
