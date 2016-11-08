PRO transform_weather_data, lores = lores

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TRANSFORMS DAILY WEATHER DATA FROM GLOBAL GRID TO
; LAND-ONLY GRID
;
; AUTHOR:
; Wolfgang Knorr, April 2008
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
dir = '/Users/wolfgang/Models/CCDAS/input/climate/'
missing = 1e20
year0 = 1979
year1 = 2005
namelist = ['cloud','precip','swdown','tmax','tmin']
units = ['fraction','mm','W/m^2','deg C','deg C']
longnames = ['cloud fraction','daily preciptiation','daily average shortwave downwelling radiation','daily maximum temperature','daily minimum temperature']
restag = 'hires'
IF KEYWORD_SET (lores) THEN restag = 'lores'
timetag = strcompress(year0,/rem)+'-'+strcompress(year1,/rem)
timetag_out = strcompress(year0,/rem)+'_'+strcompress(year1,/rem)

FOR i = 0, n_elements (namelist) - 1 DO BEGIN

; get the data
  file = dir + 'forcing.'+restag+'.'+namelist[i]+'.'+timetag+'.nc'
  ncdf_read, file = file, grid
  ncdf_read, file = file, data, variables = namelist[i], /no_struct

; construct the new grid
  land = data[*,*,0] NE missing
  l = where (reverse(land,2))
  y = n_elements(grid.lat) - 1 - l / n_elements(grid.lon)
  x = l mod n_elements(grid.lon)
  lon = grid.lon[x]
  lat = grid.lat[y]
  n = n_elements(l)
  nt = n_elements (grid.time)

; transform the format
  data_out = fltarr (nt, n)
  for j = 0, n-1 do data_out[*,j] = data[x[j],y[j],*]

; writing out the data to NetCDF
  out_file = dir+'obs_'+namelist[i]+'_daily_'+restag+'_'+timetag_out+'.nc'
  print, 'Saving to ' + out_file

  id = NCDF_CREATE(out_file, /CLOBBER)

  ; define dimensions
  dtime  = NCDF_DIMDEF(id, 'time', nt)
  dgrid  = NCDF_DIMDEF(id, 'gridnumber', n)

  vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
  NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' + $
    strcompress(year0,/rem), /CHAR
  NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR

  vgrid = NCDF_VARDEF(id, 'gridnumber', [dgrid], /FLOAT)
  NCDF_ATTPUT, id, vgrid, 'long_name', 'grid number', /CHAR

  ; define variables
  vlon = NCDF_VARDEF(id, 'lon', [dgrid], /FLOAT)
  NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
  NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

  vlat = NCDF_VARDEF(id, 'lat', [dgrid], /FLOAT)
  NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
  NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

  vxcoord = NCDF_VARDEF(id, 'xcoord', [dgrid], /LONG)
  NCDF_ATTPUT, id, vxcoord, 'long_name', $
    'integer horizontal grid coordinate', /CHAR

  vycoord = NCDF_VARDEF(id, 'ycoord', [dgrid], /LONG)
  NCDF_ATTPUT, id, vycoord, 'long_name', $
    'integer vertical grid coordinate', /CHAR

  vdata = NCDF_VARDEF(id, namelist[i], [dtime,dgrid], /FLOAT)
  NCDF_ATTPUT, id, vdata, 'long_name', longnames[i], /CHAR
  NCDF_ATTPUT, id, vdata, 'units', units[i], /CHAR

  ; define global attributes
  NCDF_ATTPUT, id, /GLOBAL, 'title', $
    'daily observed weather data for years '$
    + strcompress(year0,/rem) + ' to ' + strcompress(year1,/rem), /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'April 2008', /CHAR

  ; change to DATA mode and put data into file
  NCDF_CONTROL, id, /ENDEF

  NCDF_VARPUT, id, vlon, lon
  NCDF_VARPUT, id, vlat, lat
  NCDF_VARPUT, id, vxcoord, x
  NCDF_VARPUT, id, vycoord, y
  NCDF_VARPUT, id, vtime, indgen(nt)
  NCDF_VARPUT, id, vdata, data_out
  ; close file
  NCDF_CLOSE, id

ENDFOR

END

