PRO fapar_error_model, derr = derr, lores = lores

; generate new global FAPAR fields with explicit error
; OPTIONS:
; lores        output in lores mode
;              (default: hires)
; derr         dynamic error model, based on monthly precipitation
;              (default: simple error model, 0.1 if available)

derr_add = 0.045 ; additional error per 1mm/day monthly mean precipiation (for case 'derr')
; i.e. add error of 0.15 for every 100mm/30days precip.
;
ncdf_read,file='../input/fapar_global/fapar_hires.nc',a,/all,att=aa
fapar=a.fapar
errfapar=a.fapar*0.+0.1 ; new array for error bar in FAPAR set to 0.1
IF KEYWORD_SET (derr) THEN BEGIN
  ncdf_read,file='../input/climate/forcing.hires.monthly.1979-2008.nc',c,/all,att=ac
  t=fltarr(n_elements(a.time)) ; mapping time dimension from climate to FAPAR
  for i=0,n_elements(a.time)-1 do t[i]=where(abs(c.time-a.time[i]) lt 0.01)
  ; 
  for i=0,n_elements(t)-1 do errfapar[*,*,i]=errfapar[*,*,i]+derr_add*c.precip[*,*,t[i]]
ENDIF
errfapar(where(a.fapar lt 0.))=1e6 ; set missing values to large error (including negative FAPAR!)
lon = a.lon
lat = a.lat
;
IF KEYWORD_SET (lores) THEN BEGIN
  ncdf_read, file = '../control_bethy/adbethy_loresgrid3veg_79-00.nc', m, /all
  lsm = m.gridnum gt 0
  ni = n_elements (lon)
  nj = n_elements (lat)
  makearea, area, dlon = 360./ni, dlat = 180./nj
  lon = m.lon
  lat = m.lat
  nx = n_elements (lon)
  ny = n_elements (lat)
  nt = n_elements (a.time)
  ; transform FAPAR to lores (TM2) grid
  f2 = a.fapar ne -1
  f2 = rebin2 ((a.fapar>0)*rebin(area,ni,nj,nt)*(a.fapar ge 0.), nx, 1)
  f2 = rebin2 (f2, (ny-1)*2, 2)
  fapar = fltarr (nx, ny, nt)
  fapar[*,1:ny-2,*] = rebin (f2[*,1:(ny-2)*2,*], nx, ny-2, nt)
  fapar[*,0,*] = f2[*,0,*]
  fapar[*,ny-1,*] = f2[*,(ny-1)*2-1,*]
  a2 = rebin2 (rebin(area,ni,nj,nt)*(a.fapar ge 0.), nx, 1)
  a2 = rebin2 (a2, (ny-1)*2, 2)
  a3 = fltarr (nx, ny, nt)
  a3[*,1:ny-2,*] = rebin (a2[*,1:(ny-2)*2,*], nx, ny-2, nt)
  a3[*,0,*] = a2[*,0,*]
  a3[*,ny-1,*] = a2[*,(ny-1)*2-1,*]
  fapar = fapar / (a3 + (a3 eq 0)) * rebin (lsm, nx, ny, nt) - rebin (lsm eq 0, nx, ny, nt)
  ; transform error in FAPAR to lores (TM2) grid
  e2 = rebin2 (errfapar*rebin(area,ni,nj,nt)*(a.fapar ge 0.), nx, 1)
  e2 = rebin2 (e2, (ny-1)*2, 2)
  errfapar = fltarr (nx, ny, nt)
  errfapar[*,1:ny-2,*] = rebin (e2[*,1:(ny-2)*2,*], nx, ny-2, nt)
  errfapar[*,0,*] = e2[*,0,*]
  errfapar[*,ny-1,*] = e2[*,(ny-1)*2-1,*]
  errfapar = errfapar / (a3 + (a3 eq 0)) * rebin (lsm, nx, ny, nt) + 1e6 * rebin (lsm eq 0, nx, ny, nt)
  errfapar = errfapar * (a3 ne 0) + 1e6 * (a3 eq 0) ; cases where no available data are available
ENDIF
;
;loacdt,41
;tvf,errfapar[*,*,0],/bar,mis=1e6
;tvf,fapar[*,*,0],/bar,mis=1e6
;
; create the NetCDF output file
   outfile='../input/fapar_global/fapar_hires_serr.nc'
   IF KEYWORD_SET (lores) THEN outfile='../input/fapar_global/fapar_lores_serr.nc'
   IF KEYWORD_SET (derr) THEN outfile='../input/fapar_global/fapar_hires_derr.nc'
   IF KEYWORD_SET (derr) AND KEYWORD_SET (lores) THEN outfile='../input/fapar_global/fapar_lores_derr.nc'
   id = NCDF_CREATE(outfile, /CLOBBER)

; define dimensions
   dlon   = NCDF_DIMDEF(id, 'lon', n_elements(fapar[*,0,0]))
   dlat   = NCDF_DIMDEF(id, 'lat', n_elements(fapar[0,*,0]))
   dtime  = NCDF_DIMDEF(id, 'time', n_elements(a.time))

; description of first dimension
   vlon = NCDF_VARDEF(id, 'lon', [dlon], /FLOAT)
   NCDF_ATTPUT, id, vlon, 'long_name', aa.lon.long_name, /CHAR
   NCDF_ATTPUT, id, vlon, 'units', aa.lon.units, /CHAR

; description of second dimension
   vlat = NCDF_VARDEF(id, 'lat', [dlat], /FLOAT)
   NCDF_ATTPUT, id, vlat, 'long_name', aa.lat.long_name, /CHAR
   NCDF_ATTPUT, id, vlat, 'units', aa.lat.units, /CHAR

; description of third dimension
   vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
   NCDF_ATTPUT, id, vtime, 'long_name', aa.time.long_name, /CHAR
   NCDF_ATTPUT, id, vtime, 'units', aa.time.units, /CHAR

; define variables
   vvar1 = NCDF_VARDEF(id, 'fapar', [dlon, dlat, dtime], /FLOAT)
   NCDF_ATTPUT, id, vvar1, 'long_name', aa.fapar.long_name, /CHAR
   NCDF_ATTPUT, id, vvar1, 'units', aa.fapar.units, /CHAR
   vvar2 = NCDF_VARDEF(id, 'errfapar', [dlon, dlat, dtime], /FLOAT)
   NCDF_ATTPUT, id, vvar2, 'long_name', 'standard error', /CHAR
   NCDF_ATTPUT, id, vvar2, 'units', aa.fapar.units, /CHAR

; change to DATA mode
   NCDF_CONTROL, id, /ENDEF

; put data into file
   NCDF_VARPUT, id, vlon, lon
   NCDF_VARPUT, id, vlat, lat
   NCDF_VARPUT, id, vtime, a.time
   NCDF_VARPUT, id, vvar1, fapar
   NCDF_VARPUT, id, vvar2, errfapar

; close file
NCDF_CLOSE, id

END


