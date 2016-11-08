PRO daily_to_monthly

; converts daily input climate to monthly data

;device,decompose=0
;loadct,41

; parameters
year0_obs = 1979
year1_obs = 2008
missing = 1e6
smissing = '1e6'
;dir = '/Users/wolfgang/Models/CCDAS/input/climate/'
dir = '/Volumes/Work/ccdas_input/climate/'

; calendar stuff
dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

; input file names
tag = strcompress(year0_obs,/rem)+'-'+strcompress(year1_obs,/rem)+'.nc'
pfile = dir + 'forcing.hires.precip.'+tag
tnfile = dir + 'forcing.hires.tmin.'+tag
txfile = dir + 'forcing.hires.tmax.'+tag
rfile = dir + 'forcing.hires.swdown.'+tag
cfile = dir + 'forcing.hires.cloud.'+tag

; output file name
wfile = dir + 'forcing.hires.monthly.'+tag

; load axis information and attributes
ncdf_read, info, file = pfile, attributes = attributes
nt = n_elements(info.time)
nx = n_elements(info.lon)
ny = n_elements(info.lat)
missing_in = attributes.precip.missing_value

; generate calendar data (running month: 'm')
m = intarr(nt)
d = 0 ; day within year
mth = 0 ; month of the year
yr = year0_obs  ; calendar year
m[0] = 0 ; continuously counted month
FOR i=1,nt-1 DO BEGIN
  d = d + 1
  m[i] = m[i-1]
  dpm0 = dpm[mth]
  IF mth EQ 1 AND (yr MOD 4) EQ 0 THEN dpm0 = dpm0 + 1
  IF d EQ dpm0 THEN BEGIN
    d = 0
    IF mth EQ 11 THEN yr = yr + 1
    mth = (mth + 1) MOD 12
    m[i] = m[i] + 1
  ENDIF
ENDFOR
; check if last month in record is incomplete
IF fix(total(m EQ m[nt-1])) NE dpm0 THEN BEGIN
  ; discard last month
  lm = m[nt-1]
  m = m[where(m LT lm)]
  nt = n_elements (m)
ENDIF
nm = max(m) + 1 ; number of complete months in record

; load precipitation data
ncdf_read, precip, file = pfile, /no_struct, variables = ['precip']
precip = precip[*,*,0:nt-1]
mask = precip[*,*,0] NE missing_in

; convert to monthly mean precipitation in mm/day
mp = fltarr (nx, ny, nm)
ndays = fltarr(nm)
for i = 0, nt-1 do mp[*,*,m[i]] = mp[*,*,m[i]] + precip[*,*,i]
for i = 0, nt-1 do ndays[m[i]] = ndays[m[i]] + 1
for i = 0, nm-1 do mp[*,*,i] = mp[*,*,i] / ndays[i]
mp = mp * rebin (mask, nx, ny, nm) $
     + missing * rebin (mask EQ 0, nx, ny, nm)
a = total(temporary(precip)) ; this deletes the variable precip
print, 'done precip'

; load and convert temperature data
ncdf_read, tmax, file = txfile, /no_struct, variables = ['tmax']
tmax = tmax[*,*,0:nt-1]
mtx = fltarr (nx, ny, nm)
for i = 0, nt-1 do mtx[*,*,m[i]] = mtx[*,*,m[i]] + tmax[*,*,i]
for i = 0, nm-1 do mtx[*,*,i] = mtx[*,*,i] / ndays[i]
mtx = mtx * rebin (mask, nx, ny, nm) $
     + missing * rebin (mask EQ 0, nx, ny, nm)
a = total(temporary(tmax))
print, 'done tmax'
;
ncdf_read, tmin, file = tnfile, /no_struct, variables = ['tmin']
tmin = tmin[*,*,0:nt-1]
mtn = fltarr (nx, ny, nm)
for i = 0, nt-1 do mtn[*,*,m[i]] = mtn[*,*,m[i]] + tmin[*,*,i]
for i = 0, nm-1 do mtn[*,*,i] = mtn[*,*,i] / ndays[i]
mtn = mtn * rebin (mask, nx, ny, nm) $
     + missing * rebin (mask EQ 0, nx, ny, nm)
a = total(temporary(tmin))
print, 'done tmin'

; load and convert radiation data
ncdf_read, swdown, file = rfile, /no_struct, variables = ['swdown']
swdown = swdown[*,*,0:nt-1]
mrad = fltarr (nx, ny, nm)
for i = 0, nt-1 do mrad[*,*,m[i]] = mrad[*,*,m[i]] + swdown[*,*,i]
for i = 0, nm-1 do mrad[*,*,i] = mrad[*,*,i] / ndays[i]
mrad = mrad * rebin (mask, nx, ny, nm) $
     + missing * rebin (mask EQ 0, nx, ny, nm)
a = total(temporary(swdown))
print, 'done swdown'
;
ncdf_read, cloud, file = cfile, /no_struct, variables = ['cloud']
cloud = cloud[*,*,0:nt-1]
mcl = fltarr (nx, ny, nm)
for i = 0, nt-1 do mcl[*,*,m[i]] = mcl[*,*,m[i]] + cloud[*,*,i]
for i = 0, nm-1 do mcl[*,*,i] = mcl[*,*,i] / ndays[i]
mcl = mcl * rebin (mask, nx, ny, nm) $
     + missing * rebin (mask EQ 0, nx, ny, nm)
a = total(temporary(cloud))
print, 'done cloud'

id = NCDF_CREATE(wfile, /CLOBBER)

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
NCDF_ATTPUT, id, vtime, 'long_name', 'time', /CHAR
NCDF_ATTPUT, id, vtime, 'units', 'years', /CHAR

   ; define variables
  
vprecip = NCDF_VARDEF(id, 'precip', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vprecip, 'long_name', 'Monthly mean precipitation', /CHAR
NCDF_ATTPUT, id, vprecip, 'units', 'mm/day', /CHAR
NCDF_ATTPUT, id, vprecip, 'missing_value', smissing, /FLOAT

vtmax = NCDF_VARDEF(id, 'tmax', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vtmax, 'long_name', 'Monthly mean diurnal maximum surface air temperature', /CHAR
NCDF_ATTPUT, id, vtmax, 'units', 'degrees C', /CHAR
NCDF_ATTPUT, id, vtmax, 'missing_value', smissing, /FLOAT

vtmin = NCDF_VARDEF(id, 'tmin', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vtmin, 'long_name', 'Monthly mean diurnal minimum surface air temperature', /CHAR
NCDF_ATTPUT, id, vtmin, 'units', 'degrees C', /CHAR
NCDF_ATTPUT, id, vtmin, 'missing_value', smissing, /FLOAT

vswdown = NCDF_VARDEF(id, 'swdown', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vswdown, 'long_name', 'Monthly mean shortwave downwelling radiation', /CHAR
NCDF_ATTPUT, id, vswdown, 'units', 'W / m^2', /CHAR
NCDF_ATTPUT, id, vswdown, 'missing_value', smissing, /FLOAT

vcloud = NCDF_VARDEF(id, 'cloud', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vcloud, 'long_name', 'Monthly mean total cloud fraction', /CHAR
NCDF_ATTPUT, id, vcloud, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vcloud, 'missing_value', smissing, /FLOAT

   ; define global attributes
NCDF_ATTPUT, id, /GLOBAL, 'title', 'Monthly averaged gridded station data', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR

   ; change to DATA mode and put data into file
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, vlon, info.lon
NCDF_VARPUT, id, vlat, info.lat
NCDF_VARPUT, id, vtime, (indgen(nm)+0.5)/12.+year0_obs
NCDF_VARPUT, id, vprecip, mp
NCDF_VARPUT, id, vtmax, mtx
NCDF_VARPUT, id, vtmin, mtn
NCDF_VARPUT, id, vswdown, mrad
NCDF_VARPUT, id, vcloud, mcl

   ; close file
NCDF_CLOSE, id

END





