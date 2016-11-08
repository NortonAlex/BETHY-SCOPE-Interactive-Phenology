PRO extract_sites, sites = sites, plot = plot

; parameters
year0_obs = 1979
year1_obs = 2008
missing = 1e6
smissing = '1e6'
dir = '../input/climate/'
odir = '../input/eddy_sites/'

; calendar stuff
dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month_name = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

; input file names
tag = strcompress(year0_obs,/rem)+'-'+strcompress(year1_obs,/rem)
pfile = dir + 'forcing.hires.precip.'+tag+'.nc'
tnfile = dir + 'forcing.hires.tmin.'+tag+'.nc'
txfile = dir + 'forcing.hires.tmax.'+tag+'.nc'
rfile = dir + 'forcing.hires.swdown.'+tag+'.nc'
sfile = odir + 'forcing_'+tag+'_'
;
grid_dir = '../control_bethy/'
grid_file = grid_dir + 'bethy_grid3veg_v2.nc'
;
sites_dir = '../control_bethy/'
sites_file = sites_dir + 'site_specs_all.dat'
IF KEYWORD_SET (sites) THEN sites_file = sites_dir + sites
;
; get co-ordinates and names of sites
max_sites = 1000
site_name = strarr (max_sites)
lon_site = fltarr (max_sites)
lat_site = fltarr (max_sites)
openr, u, sites_file, /get_lun
line = ' '
readf, u, line ; skip the first line
FOR j = 0, max_sites - 1 DO BEGIN
  IF eof(u) THEN GOTO, EXIT_LOOP
  readf, u, line
  reads, line, xlat, xlon
  k = strpos (line, '!')
  site_name[j] = strcompress(strmid(line,k+1),/rem)
  lon_site[j] = xlon
  lat_site[j] = xlat
ENDFOR
EXIT_LOOP:
free_lun, u
;jsite = jsite[0:j-1]
site_name = site_name[0:j-1]
lon_site = lon_site[0:j-1]
lat_site = lat_site[0:j-1]
nsite = j

; map onto main grid
ncdf_read, file=grid_file, grid, variables=['lat', 'lon', 'gridnum']
l = where (grid.gridnum gt 0) ; land grid
x = l mod n_elements(grid.lon)
y = l / n_elements(grid.lon)
jsite = lonarr (nsite)
FOR j = 0, nsite - 1 DO BEGIN
  l = sort (abs(grid.lon[x]-lon_site[j])+abs(grid.lat[y]-lat_site[j]))
  jsite[j] = l[0]
ENDFOR
; checking:
IF KEYWORD_SET (plot) THEN BEGIN
  plot,grid.lon[x],grid.lat[y],psym=3,xrange=[-180,180],/xst,yrange=[-90,90],/yst
  oplot,grid.lon[x[jsite]],grid.lat[y[jsite]],psym=2
ENDIF
FOR j = 0, nsite - 1 DO $
  PRINT, site_name[j], lat_site[j], grid.lat[y[jsite[j]]], $
  lon_site[j], grid.lon[x[jsite[j]]]

; load the precipitation data
ncdf_read, precip, file = pfile, /no_struct, variables = ['precip']
nt = n_elements(precip[0,0,*])
; re-map onto site locations only
p = fltarr (nt, nsite)
FOR j = 0, nsite-1 DO p[*,j] = precip[x[jsite[j]],y[jsite[j]],*]
a = total(temporary(precip)) ; this deletes the variable precip

; load and map the temperature data
ncdf_read, tmax, file = txfile, /no_struct, variables = ['tmax']
tx = fltarr (nt, nsite) & IF n_elements(tmax[0,0,*]) NE nt THEN STOP
FOR j = 0, nsite-1 DO tx[*,j] = tmax[x[jsite[j]],y[jsite[j]],*]
a = total(temporary(tmax))
ncdf_read, tmin, file = tnfile, /no_struct, variables = ['tmin']
tn = fltarr (nt, nsite) & IF n_elements(tmin[0,0,*]) NE nt THEN STOP
FOR j = 0, nsite-1 DO tn[*,j] = tmin[x[jsite[j]],y[jsite[j]],*]
a = total(temporary(tmin))

; load the incoming shortwave radiation data
ncdf_read, swdown, file = rfile, /no_struct, variables = ['swdown']
rad = fltarr (nt, nsite) & IF n_elements(swdown[0,0,*]) NE nt THEN STOP
FOR j = 0, nsite - 1 DO rad[*,j] = swdown[x[jsite[j]],y[jsite[j]],*]
a = total(temporary(swdown))

  ; save in one file per site
  FOR j = 0, nsite-1 DO BEGIN

    id = NCDF_CREATE(sfile+site_name[j]+'.nc', /CLOBBER)

    ; define dimensions
    dtime  = NCDF_DIMDEF(id, 'time', nt)

    vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
    NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' + $
      strcompress(year0_obs,/rem), /CHAR
    NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR

    ; define variables
    vlon = NCDF_VARDEF(id, 'lon', /FLOAT)
    NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
    NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

    vlat = NCDF_VARDEF(id, 'lat', /FLOAT)
    NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
    NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

    vprecip = NCDF_VARDEF(id, 'precip', [dtime], /FLOAT)
    NCDF_ATTPUT, id, vprecip, 'long_name', 'Daily precipitation', /CHAR
    NCDF_ATTPUT, id, vprecip, 'units', 'mm', /CHAR
    NCDF_ATTPUT, id, vprecip, 'missing_value', smissing, /FLOAT

    vtmin = NCDF_VARDEF(id, 'tmin', [dtime], /FLOAT)
    NCDF_ATTPUT, id, vtmin, 'long_name', 'Daily minimum temperatue', /CHAR
    NCDF_ATTPUT, id, vtmin, 'units', 'degrees Celsius', /CHAR
    NCDF_ATTPUT, id, vtmin, 'missing_value', smissing, /FLOAT

    vtmax = NCDF_VARDEF(id, 'tmax', [dtime], /FLOAT)
    NCDF_ATTPUT, id, vtmax, 'long_name', 'Daily minimum temperatue', /CHAR
    NCDF_ATTPUT, id, vtmax, 'units', 'degrees Celsius', /CHAR
    NCDF_ATTPUT, id, vtmax, 'missing_value', smissing, /FLOAT

    vswdown = NCDF_VARDEF(id, 'swdown', [dtime], /FLOAT)
    NCDF_ATTPUT, id, vswdown, 'long_name', $
      'Daily mean incoming solar radation', /CHAR
    NCDF_ATTPUT, id, vswdown, 'units', 'W / m^2', /CHAR
    NCDF_ATTPUT, id, vswdown, 'missing_value', smissing, /FLOAT

    ; define global attributes
    NCDF_ATTPUT, id, /GLOBAL, 'title', $
      'daily generated climate data for year '$
      + strcompress(year0_obs,/rem) + ' to ' + strcompress(year1_obs,/rem), /CHAR
    NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
    NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
    NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'February 2010', /CHAR

    ; change to DATA mode and put data into file
    NCDF_CONTROL, id, /ENDEF

    NCDF_VARPUT, id, vlon, lon_site[j]
    NCDF_VARPUT, id, vlat, lat_site[j]
    NCDF_VARPUT, id, vtime, indgen(nt)
    NCDF_VARPUT, id, vprecip, p[*,j]
    NCDF_VARPUT, id, vtmin, tn[*,j]
    NCDF_VARPUT, id, vtmax, tx[*,j]
    NCDF_VARPUT, id, vswdown, rad[*,j]
    ; close file
    NCDF_CLOSE, id

  ENDFOR

END



