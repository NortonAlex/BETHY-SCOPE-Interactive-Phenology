PRO adjust_pft, lores = lores, ratio = ratio, nowindow=nowindow
;
; Adjust PFT fractions in CCDAS grid file such that prior FAPAR average matches observations
; or simply calculate ratio defined as temporal mean FAPAR of observations devided by
; temporal mean of modelled FAPAR and store in grid file as additional data field
; Expects result from post run in file 'fapar_global.dat' in directory '../output/'
; in correct resolution for years 1999 to 2005, and
; observed FAPAR from May 2002 to December 2003
; -> adapt directory of satellite observed FAPAR
; -> adjust colour table number if necessary
;
; OPTIONS:
; lores         for low resolution (default: high resolution)
; ratio         save ratio instead of adjusting PFT fractions
;
fapar_model = '../output/fapar_global.dat'
; take only May 2002 to December 2005 (last 44 months) for 7-year run:
nyears = 7
last_month_model = nyears*12-1
first_month_model = last_month_model-44+1
fapar_obs = '~/Data/MERIS/MERIS_global/fapar_hires.nc'
IF KEYWORD_SET (lores) THEN fapar_obs = '~/Data/MERIS/MERIS_global/fapar_lores.nc'
oldCCDASfile = '../control_bethy/bethy_grid3veg_v2.nc'
IF KEYWORD_SET (lores) THEN $
  oldCCDASfile = '../control_bethy/bethy_loresgrid3veg_v2.nc'
newCCDAS_outfile = '../control_bethy/bethy_grid3veg_pftcorr.nc'
IF KEYWORD_SET (lores) THEN newCCDAS_outfile = $
  '../control_bethy/bethy_loresgrid3veg_pftcorr.nc'
IF KEYWORD_SET (ratio) THEN newCCDAS_outfile = $
  '../control_bethy/bethy_grid3veg_faparratio.nc'
IF KEYWORD_SET (ratio) AND KEYWORD_SET (lores) THEN newCCDAS_outfile = $
  '../control_bethy/bethy_loresgrid3veg_faparratio.nc'
;
dlonlat = [2., 2.]
IF KEYWORD_SET (lores) THEN dlonlat = [10.,7.5]
; 
; model FAPAR lores for 5 years (1999 to 2003):
loadm,'../output/fpar_global.dat',fpar,year=-1,gr=dlonlat,/center,oc=-1
fpar=fpar[*,*,first_month_model:last_month_model]
; MERIS FAPAR lores:
ncdf_read,file=fapar_obs,d,/all
;
; pixel-by-pixel bias:
corr=total(d.fapar,3)/total(fpar,3)<1
;
; visual check
n=last_month_model-first_month_model+1
IF NOT KEYWORD_SET (nowindow) THEN BEGIN
loadct,42
window,0,xs=720,ys=460
tvf,total(d.fapar,3)/n,/bar,mis=-1,min=0,max=1,title='mean satellite FAPAR'
window,1,xs=720,ys=460
tvf,total(fpar,3)/n,/bar,mis=-1,min=0,max=1,title='mean model FAPAR'
window,2,xs=720,ys=460
tvf,corr,/bar,title='correction factor'
ENDIF
;
; load PFT fractions
ncdf_read,file=oldCCDASfile,b,/all
nx = n_elements(b.lon)
ny = n_elements(b.lat)
npft = n_elements(b.type[0,0,*])
IF NOT KEYWORD_SET (ratio) THEN $
  b.frac=b.frac*rebin(corr,nx,ny,npft)
;tvf,total(b.frac,3),miss=-9999.*npft,/bar

; create new netcdf file:
id = NCDF_CREATE(newCCDAS_outfile, /CLOBBER)

  ; define dimensions
dlon   = NCDF_DIMDEF(id, 'lon', n_elements(b.lon))
dlat   = NCDF_DIMDEF(id, 'lat', n_elements(b.lat))
dvnum  = NCDF_DIMDEF(id, 'vnum', 3)
dtime  = NCDF_DIMDEF(id, 'time', /UNLIMITED)


vlon = NCDF_VARDEF(id, 'lon', [dlon], /FLOAT)
NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

vlat = NCDF_VARDEF(id, 'lat', [dlat], /FLOAT)
NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

vvnum =  NCDF_VARDEF(id, 'vnum', [dvnum], /LONG)
NCDF_ATTPUT, id, vvnum, 'long_name', 'vegetation number', /CHAR
NCDF_ATTPUT, id, vvnum, 'units', '', /CHAR

vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
NCDF_ATTPUT, id, vtime, 'long_name', 'time', /CHAR
NCDF_ATTPUT, id, vtime, 'units', 'month', /CHAR

   ; define variables
  
vgridnum = NCDF_VARDEF(id, 'gridnum', [dlon, dlat], /LONG)
NCDF_ATTPUT, id, vgridnum, 'long_name', 'land grid cell number starting at 180E/90N', /CHAR
NCDF_ATTPUT, id, vgridnum, 'units', '', /CHAR
NCDF_ATTPUT, id, vgridnum, 'missing_value', '-9999', /LONG

velev = NCDF_VARDEF(id, 'elev', [dlon, dlat], /LONG)
NCDF_ATTPUT, id, velev, 'long_name', 'Mean grid cell elevation', /CHAR
NCDF_ATTPUT, id, velev, 'units', 'm', /CHAR
NCDF_ATTPUT, id, velev, 'missing_value', '-9999', /LONG

vtype = NCDF_VARDEF(id, 'type', [dlon, dlat, dvnum], /LONG)
NCDF_ATTPUT, id, vtype, 'long_name', 'Vegetation type (0-13)', /CHAR
NCDF_ATTPUT, id, vtype, 'units', '', /CHAR
NCDF_ATTPUT, id, vtype, 'missing_value', '-9999', /LONG

vfrac = NCDF_VARDEF(id, 'frac', [dlon, dlat, dvnum], /FLOAT)
NCDF_ATTPUT, id, vfrac, 'long_name', 'Vegetation type fraction', /CHAR
NCDF_ATTPUT, id, vfrac, 'units', '', /CHAR
NCDF_ATTPUT, id, vfrac, 'missing_value', '-9999', /FLOAT

vwmax = NCDF_VARDEF(id, 'wmax',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vwmax, 'long_name', 'Maximum plant-available soil moisture', /CHAR
NCDF_ATTPUT, id, vwmax, 'units', 'mm', /CHAR
NCDF_ATTPUT, id, vwmax, 'missing_value', '-9999', /FLOAT

vrhosw = NCDF_VARDEF(id, 'rhosw',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vrhosw, 'long_name', 'Wet-soil albedo', /CHAR
NCDF_ATTPUT, id, vrhosw, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vrhosw, 'missing_value', '-9999', /FLOAT

vrhosd = NCDF_VARDEF(id, 'rhosd',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vrhosd, 'long_name', 'Dry-soil albedo', /CHAR
NCDF_ATTPUT, id, vrhosd, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vrhosd, 'missing_value', '-9999', /FLOAT

vdes = NCDF_VARDEF(id, 'des',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vdes, 'long_name', 'Soil water desortptivity', /CHAR
NCDF_ATTPUT, id, vdes, 'units', 'mm/sqrt(day)', /CHAR
NCDF_ATTPUT, id, vdes, 'missing_value', '-9999', /FLOAT

vevap1 = NCDF_VARDEF(id, 'evap1',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vevap1, 'long_name', 'Phase 1 soil evaporation', /CHAR
NCDF_ATTPUT, id, vevap1, 'units', 'mm', /CHAR
NCDF_ATTPUT, id, vevap1, 'missing_value', '-9999', /FLOAT

vmpot = NCDF_VARDEF(id, 'mpot', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vmpot, 'long_name', 'Mean potential shortwave radiation', /CHAR
NCDF_ATTPUT, id, vmpot, 'units', 'W/m2', /CHAR
NCDF_ATTPUT, id, vmpot, 'missing_value', '-9999', /FLOAT

IF KEYWORD_SET (ratio) THEN BEGIN
vratio = NCDF_VARDEF(id, 'ratio', [dlon, dlat], /FLOAT)
NCDF_ATTPUT, id, vratio, 'long_name', 'mean observed / mean prior model FAPAR', /CHAR
NCDF_ATTPUT, id, vratio, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vratio, 'missing_value', '-9999', /FLOAT
ENDIF

   ; define global attributes
NCDF_ATTPUT, id, /GLOBAL, 'title', 'Boundary conditions for CCDAS-BETHY', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'University of Bristol', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'date', 'Version with bias correced PFT fractions, creation date 31 March 2010', /CHAR

   ; change to DATA mode and put data into file
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, vlon, b.lon
NCDF_VARPUT, id, vlat, b.lat
NCDF_VARPUT, id, vvnum, b.vnum
NCDF_VARPUT, id, vtime, b.time
NCDF_VARPUT, id, vgridnum, b.gridnum
NCDF_VARPUT, id, velev, b.elev
NCDF_VARPUT, id, vtype, b.type
NCDF_VARPUT, id, vfrac, b.frac
NCDF_VARPUT, id, vwmax, b.wmax
NCDF_VARPUT, id, vrhosw, b.rhosw
NCDF_VARPUT, id, vrhosd, b.rhosd
NCDF_VARPUT, id, vdes, b.des
NCDF_VARPUT, id, vevap1, b.evap1
NCDF_VARPUT, id, vmpot, b.mpot
IF KEYWORD_SET (ratio) THEN $
  NCDF_VARPUT, id, vratio, corr*(b.gridnum ge 0)-9999*(b.gridnum lt 0)

   ; close file
NCDF_CLOSE, id
print, 'Written '+newCCDAS_outfile

;STOP

; TEST:
IF NOT KEYWORD_SET (nowindow) THEN BEGIN
IF KEYWORD_SET (ratio) THEN BEGIN
ncdf_read,file=newCCDAS_outfile,b2,/all,attr=a
window,3,xs=720,ys=460
tvf,b2.ratio,miss=a.frac.missing_value,/bar,title='BIAS RATIO IN FILE'
ENDIF ELSE BEGIN
ncdf_read,file=newCCDAS_outfile,b2,/all,attr=a
window,3,xs=720,ys=460
tvf,total(b2.frac,3),/bar,miss=a.frac.missing_value*3,title='NEW TOTAL PFT FRACTION'
ncdf_read,file=oldCCDASfile,b2,/all,attr=a
window,4,xs=720,ys=460
tvf,total(b2.frac,3),/bar,miss=a.frac.missing_value*3,title='OLD TOTAL PFT FRACTION'
ENDELSE
ENDIF

END



