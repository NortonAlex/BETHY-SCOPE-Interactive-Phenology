PRO daily_climate_parameters, sites = sites, lores = lores, screenplot = screenplot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (1) INITIATALIZATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, 'Initiatlizing...'

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
month_name = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

; input file names
tag = strcompress(year0_obs,/rem)+'-'+strcompress(year1_obs,/rem)+'.nc'
pfile = dir + 'forcing.hires.precip.'+tag
tnfile = dir + 'forcing.hires.tmin.'+tag
txfile = dir + 'forcing.hires.tmax.'+tag
rfile = dir + 'forcing.hires.swdown.'+tag
cfile = dir + 'forcing.hires.cloud.'+tag
;
;grid_dir = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/'
grid_dir = '/Volumes/Work/wolfgang/CCDAS/ppm/control_bethy/'
grid_file = grid_dir + 'bethy_grid3veg_v2.nc'
grid_file_lores = grid_dir + 'bethy_loresgrid3veg_v2.nc'
;
;sites_dir = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/'
sites_dir = '/Volumes/Work/wolfgang/CCDAS/ppm/control_bethy/'
sites_file = sites_dir + 'site_specs_7s.dat'

; load elevation and grid information
ncdf_read, file = grid_file, variables = ['gridnum', 'elev'], grid
land = grid.gridnum ne -9999.
; the grid cells of the basic grid are numbered S to N and E to W
l = where (land)
y = l / n_elements(grid.lon)
x = l mod n_elements(grid.lon)

IF KEYWORD_SET (sites) THEN BEGIN
  ; get co-ordinates of sites
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
ENDIF ELSE IF KEYWORD_SET (lores) THEN BEGIN
  ncdf_read, file = grid_file_lores, variables = ['gridnum'], grid0
  land0 = grid0.gridnum ne -9999.
  l = where (reverse(land0,2))
  y0 = n_elements(grid0.lat) - 1 - l / n_elements(grid0.lon)
  x0 = l mod n_elements(grid0.lon)
  lon_site = grid0.lon[x0]
  lat_site = grid0.lat[y0]
  nsite = n_elements(l)
ENDIF ELSE BEGIN
  l = where (reverse(land,2))
  y0 = n_elements(grid.lat) - 1 - l / n_elements(grid.lon)
  x0 = l mod n_elements(grid.lon)
  lon_site = grid.lon[x0]
  lat_site = grid.lat[y0]
  nsite = n_elements(l)
ENDELSE

; map onto main grid
jsite = lonarr (nsite)
FOR j = 0, nsite - 1 DO BEGIN
  l = sort (abs(grid.lon[x]-lon_site[j])+abs(grid.lat[y]-lat_site[j]))
  jsite[j] = l[0]
ENDFOR
; checking:
IF KEYWORD_SET (screenplot) THEN BEGIN
  plot,grid.lon[x],grid.lat[y],psym=3,xrange=[-180,180],/xst,yrange=[-90,90],/yst
  oplot,grid.lon[x[jsite]],grid.lat[y[jsite]],psym=2
ENDIF

; the monthly climate statistics
mp = fltarr(12,nsite) ; precipitation...
mw = fltarr(12,nsite)
mww = fltarr(12,nsite)
mwd = fltarr(12,nsite)
mpkwb = fltarr(12,nsite) ; shape parameter of Weibull distribution
mtmw = fltarr(12,nsite) ; temperature...
mtsw = fltarr(12,nsite)
mtrmw = fltarr(12,nsite)
mtrsw = fltarr(12,nsite)
mtmd = fltarr(12,nsite)
mtsd = fltarr(12,nsite)
mtrmd = fltarr(12,nsite)
mtrsd = fltarr(12,nsite)
mrmw = fltarr(12,nsite) ; radation...
mrsw = fltarr(12,nsite)
mrmd = fltarr(12,nsite)
mrsd = fltarr(12,nsite)
; for quadratic fit, actual to pot. radiation ratio as a function of cloudiness
aratio = fltarr(12,3,nsite)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (2) STATISTICS FOR PRECIPITATION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, 'Processing precipitation data...'

; load the precipitation data
ncdf_read, precip, file = pfile, /no_struct, variables = ['precip']
nt = n_elements(precip[0,0,*])
precip_min = min (precip[where(precip GT 0. AND precip NE missing)])
print, 'minimum precipitation in data set: ', precip_min

; generate calendar data for observations (day, month, year)
yr = intarr(nt) & m = intarr(nt) & d = intarr(nt) & jday = intarr(nt)
d[0] = 0 & m[0] = 0 & yr[0] = 0 & jday[0] = 1
FOR i=1,nt-1 DO BEGIN
  d[i] = d[i-1] + 1
  jday[i] = jday[i-1] + 1
  m[i] = m[i-1]
  yr[i] = yr[i-1]
  dpm0 = dpm[m[i]]
  IF m[i] EQ 1 AND ((yr[i]+year0_obs) MOD 4) EQ 0 THEN dpm0 = dpm0 + 1
  IF d[i] EQ dpm0 THEN BEGIN
    d[i] = 0
    IF m[i] EQ 11 THEN BEGIN
      yr[i] = yr[i] + 1
      jday[i] = 1
    ENDIF
    m[i] = (m[i] + 1) MOD 12
  ENDIF
ENDFOR
nyr = year1_obs - year0_obs + 1
IF max(yr)+1 NE nyr THEN STOP
l = where (m+1 EQ 12 AND d+1 EQ 31)
l = l[n_elements(l)-1] ; last 31 December in record
print, 'days in record: ', nt
print, 'days until including end of last complete year: ', l+1
; if last year is incomplete, discard it
; (conforming with 'daily_to_monthly.pro')
IF nt NE l+1 THEN BEGIN
  nt = l+1
  year1_obs = year1_obs - 1
ENDIF

; re-map onto site locations only
p = fltarr (nt, nsite)
FOR j = 0, nsite-1 DO p[*,j] = precip[x[jsite[j]],y[jsite[j]],0:nt-1]
a = total(temporary(precip)) ; this deletes the variable precip
IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  set_plot,'ps'
  device,file='precip_sites.ps'
  for j=0,nsite-1 do plot,p[*,j],title='precip '+site_name[j]
  device,/close
  set_plot,'x'
ENDIF

; compute the statistics for the rainfall at the sites
ndays = intarr (12) & for i = 0, 11 do ndays[i] = total (m EQ i)
; make look-up table for inversion of Gamma functions
;nwb = 1000000
;kmin = 0.3
;kwb = findgen(nwb+1)/nwb * (1.-kmin) + kmin
;grwb = gamma(1.+1./kwb)^2 / gamma(1.+2./kwb)
;grmin = grwb[0] & grmax = grwb[nwb]
;;plot,kwb,grwb
;ntwb = 100000 ; size of the table
;grt = findgen(ntwb+1)/ntwb*(grmax-grmin)+grmin
;kt = fltarr(ntwb+1)
;kt[0] = kmin & kt[ntwb] = 1.0
;j = 1L  ; make sure these are long integers!
;FOR i = 1L, ntwb-1 DO BEGIN
;  REPEAT j=j+1 UNTIL grwb[j] GE grt[i] ; always search forward (steady function)
;  kt[i] = kwb[j]
;  ; print,i,grt[i],j,kwb[j],grwb[j]
;  j = j - 1
;ENDFOR
FOR j = 0, nsite - 1 DO BEGIN
  ; number of wet days following wet day per month
  for i=1,nt-1 do mww[m[i],j] = mww[m[i],j] $
                                    + (p[i,j] GT 0. AND p[i-1,j] GT 0.)
  ; number of wet days following dry day per month
  for i=1,nt-1 do mwd[m[i],j] = mwd[m[i],j] $
                                  + (p[i,j] GT 0. AND p[i-1,j] EQ 0.)
  FOR i=0,11 DO BEGIN
    l0 = where (p[*,j] GT 0. AND m EQ i)
    n0 = n_elements (l0)
    ; monthly mean precip [mm/day]
    IF l0[0] NE -1 THEN mp[i,j] = total (p[l0,j])/ndays[i] $
      ELSE mp[i,j] = 0.
    ; number of wet days per month
    mw[i,j] = n0
    ; moments estimation for Weibull distribution
    IF l0[0] EQ -1 THEN BEGIN
      mpkwb[i,j] = missing
    ENDIF ELSE IF n0 LT 6 THEN BEGIN
      mpkwb[i,j] = 0.8  ; default Weibull shape parameter
    ENDIF ELSE BEGIN
      ; estimate 'k' from ratio of percentiles
      prc1 = 1./3.
      ;prc2 = 9./10.
      ;prc2 = 1. - 10./n0
      ;if n0 LT 100 then prc2 = 2./3.
      prc2 = 2./3.
      s = sort (p[l0,j])
      p1 = p[l0[s[prc1*(n0-1)]],j]
      IF prc1*(n0-1) MOD 1. GT 0. THEN $
        p1 = (p1 + p[l0[s[prc1*(n0-1)+1]],j])/2.
      p2 = p[l0[s[prc2*(n0-1)]],j]
      IF prc2*(n0-1) MOD 1. GT 0. THEN $
        p2 = (p2 + p[l0[s[prc2*(n0-1)+1]],j])/2.
      mpkwb[i,j] = alog(alog(1./(1.-prc2))/alog(1./(1.-prc1)))/alog(p2/p1)
      result = alog(alog(1./(1.-prc2))/alog(1./(1.-prc1)))/alog(p2/p1)
      IF finite(result) THEN mpkwb[i,j] = result ELSE mpkwb[i,j] = 0.8
      if keyword_set(sites) then $
        print, site_name[j], j, i, n0, p1, p2, mpkwb[i,j]
      ; estimate 'k' from ratio of moments
      ; mom1 = total(p[l0,j]-precip_min)/n0     ; first moment
      ;mom2 = total((p[l0,j]-precip_min)^2)/n0 ; second moment
      ;idx = (fix((mom1^2/mom2-grmin)/(grmax-grmin)*ntwb)+1) > 0 < ntwb
      ;mpkwb[i,j] = kt[idx]
    ENDELSE
  ENDFOR
  IF j mod 100 EQ 0 THEN print, 'Precipitation, land point ', j, '   of ', nsite
ENDFOR
; divide by number of wet days
for i=0,11 do mww[i,*] = mww[i,*] / (mw[i,*] + (mw[i,*] EQ 0))
; divide by number of dry days
for i=0,11 do mwd[i,*] = mwd[i,*] / (ndays[i] - mw[i,*] + $
                           (mw[i,*] EQ ndays[i]))
for i=0,11 do mw[i,*] = mw[i,*] / ndays[i]

IF KEYWORD_SET (sites) THEN BEGIN
  ; plot the rainfall distribution with the fit (to file)
  set_plot,'ps'
  device,file='raindistribution_sites_weibull.ps'
  FOR j = 0, nsite-1 DO BEGIN
    FOR i = 0, 11 DO BEGIN
      ; histograms
      l0 = where (p[*,j] GT 0. AND m EQ i)
      nbins = n_elements(l0)/20 > 3
      h = histogram(p[l0,j],omin=omin,omax=omax,nbins=nbins)
      xx = (findgen(nbins) + 0.5) / nbins * $
          (omax-omin) + omin
      plot,xx,h,/ylog,yrange=[0.1,1e3],xrange=[omin,omax],/xst,psym=2,$
        title=month_name[i]+' rainfall at '+site_name[j]+' '+$
          strcompress(year0_obs,/rem)+' to '+strcompress(year1_obs,/rem)$
          +'  k='+strcompress(mpkwb[i,j],/rem), $
        ytitle='Number of days (line: Weibull model)',$
        xtitle='rainfall [mm/day]'
      lwb = mp[i,j] / mw[i,j] / gamma(1.+1./mpkwb[i,j])
      wb = mpkwb[i,j]/lwb*((xx-precip_min)/lwb)^(mpkwb[i,j]-1)*$
           exp(-((xx-precip_min)/lwb)^mpkwb[i,j])
      wb = wb / total(wb) * total(h)
      oplot,xx,wb
    ENDFOR
  ENDFOR
  device,/close
  set_plot,'x'
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (3) STATISTICS FOR TEMPERATURE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, 'Processing temperature data...'

; convert precipitation data to byte (rain / no-rain) to save space
p = p GT 0.

; load the temperature data
ncdf_read, tmax, file = txfile, /no_struct, variables = ['tmax']
tx = fltarr (nt, nsite)
FOR j = 0, nsite-1 DO tx[*,j] = tmax[x[jsite[j]],y[jsite[j]],0:nt-1]
a = total(temporary(tmax))
ncdf_read, tmin, file = tnfile, /no_struct, variables = ['tmin']
tn = fltarr (nt, nsite)
FOR j = 0, nsite-1 DO tn[*,j] = tmin[x[jsite[j]],y[jsite[j]],0:nt-1]
a = total(temporary(tmin))
IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  set_plot,'ps'
  device,file='temperature_sites.ps',/color
  for j=0,nsite-1 do begin & $
    plot,tx[*,j],title='Temperature '+site_name[j],color=1 & $
    oplot,tn[*,j],color=128 & endfor
  device,/close
  set_plot,'x'
ENDIF

; compute the statistics for the temperature mean and range
FOR j = 0, nsite - 1 DO BEGIN

  tm = reform(tx[*,j]+tn[*,j]) / 2.
  tr = reform(tx[*,j]-tn[*,j])

  FOR i=0,11 DO BEGIN

    l0 = where (p[*,j] GT 0 AND m EQ i)
    IF l0[0] EQ -1 THEN BEGIN
      m1 = [missing, missing^2]
      m2 = [missing, missing^2]
    ENDIF ELSE IF n_elements(l0) EQ 1 THEN BEGIN
      m1 = [tm[l0], missing^2]
      m2 = [tr[l0], missing^2]
    ENDIF ELSE BEGIN
      m1 = moment (tm[l0])
      m2 = moment (tr[l0])
    ENDELSE

    ; daily mean temperature, mean and std. dev. on rain days
    mtmw[i,j] = m1[0]
    mtsw[i,j] = sqrt(m1[1])

    ; daily temperature range, mean and std. dev. on rain days
    mtrmw[i,j] = m2[0]
    mtrsw[i,j] = sqrt(m2[1])

    l0 = where (p[*,j] EQ 0 AND m EQ i)
    IF l0[0] EQ -1 THEN BEGIN
      m1 = [missing, missing^2]
      m2 = [missing, missing^2]
    ENDIF ELSE IF n_elements(l0) EQ 1 THEN BEGIN
      m1 = [tm[l0], missing^2]
      m2 = [tr[l0], missing^2]
    ENDIF ELSE BEGIN
      m1 = moment (tm[l0])
      m2 = moment (tr[l0])
    ENDELSE

    ; daily mean temperature, mean and std. dev. on dry days
    mtmd[i,j] = m1[0]
    mtsd[i,j] = sqrt(m2[1])

    ; daily temperature range, mean and std. dev. on dry days
    mtrmd[i,j] = m2[0]
    mtrsd[i,j] = sqrt(m2[1])

    ; if missing value in rain or dry days, use from the other
    IF mtmw[i,j] EQ missing THEN mtmw[i,j] = mtmd[i,j]
    IF mtsw[i,j] EQ missing THEN mtsw[i,j] = mtsd[i,j]
    IF mtrmw[i,j] EQ missing THEN mtrmw[i,j] = mtrmd[i,j]
    IF mtrsw[i,j] EQ missing THEN mtrsw[i,j] = mtrsd[i,j]
    IF mtmd[i,j] EQ missing THEN mtmd[i,j] = mtmw[i,j]
    IF mtsd[i,j] EQ missing THEN mtsd[i,j] = mtsw[i,j]
    IF mtrmd[i,j] EQ missing THEN mtrmd[i,j] = mrtmw[i,j]
    IF mtrsd[i,j] EQ missing THEN mtrsd[i,j] = mtrsw[i,j]

  ENDFOR
  IF j mod 100 EQ 0 THEN print, 'Temperature, land point ', j, '   of ', nsite
ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (4) STATISTICS FOR RADIATION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, 'Processing radiation data...'

; compute standard atm. pressure for potential radiation
; delete temperature variable to free space
laps = 0.006
p0 = 1.01325e5
pair = fltarr (nt, nsite)
for j = 0, nsite-1 do pair[*,j] = $
    p0*(1./(1.+grid.elev[x[jsite[j]],y[jsite[j]]]*laps/$
    ((tn[*,j]+tx[*,j])/2.+273.15)))
a = total(temporary(tn))
a = total(temporary(tx))

; parameters for potential solar incoming radiation
a0 = [1.00011, 0.34221E-1, 0.128E-2, 0.719E-3, 0.77E-4]
s0 = 1360.    ;  total radiation
b0 = 0.185
d0 = 0.4

ndays0 = 365.


; compute potential radiation
; delete pressure variable to free space
rpot = fltarr (nt, nsite)
FOR i = 0, nt - 1 DO BEGIN

  alpha1 = 2. * !pi * (jday[i] - 1) / ndays0
  alpha2 = alpha1 * 2.0

  delta= -23.4 * COS (2. * !pi * (jday[i] + 10) / ndays0 )

  spd = SIN (grid.lat[y[jsite]] * !pi / 180.) * SIN (delta * !pi / 180.)
  cpd = COS (grid.lat[y[jsite]] * !pi / 180.) * COS (delta * !pi / 180.)
  dbodsq = a0(0) + a0(1) * COS (alpha1) + a0(2) * SIN (alpha1) + $
           a0(3) * COS (alpha2) + a0(4) * SIN (alpha2)

  FOR hour = 0, 23 DO BEGIN
    h = hour * !pi / 12.
    coszen =   (spd - cpd * COS (h)) > 1e-12
    rtop = s0 * dbodsq * coszen
    tdir = EXP (-b0 / coszen * pair[i,*] / p0)

    rpot[i,*] = rpot[i,*] + rtop * (d0 + (1. - d0) * tdir)
  ENDFOR
  IF i mod 1000 EQ 0 THEN print, 'Potential radiation, day ', i, '    of ', nt

ENDFOR
rpot = rpot / 24.
a = total(temporary(pair))
IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  set_plot,'ps'
  device,file='potswdown_sites.ps'
  for j=0,nsite-1 do plot,rpot[*,j],title='Potential rad. '+site_name[j],color=1
  device,/close
  set_plot,'x'
ENDIF

; load the incoming shortwave radiation data
ncdf_read, swdown, file = rfile, /no_struct, variables = ['swdown']
rad = fltarr (nt, nsite)
FOR j = 0, nsite - 1 DO rad[*,j] = swdown[x[jsite[j]],y[jsite[j]],0:nt-1]
a = total(temporary(swdown))
IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  set_plot,'ps'
  device,file='swdown_sites.ps'
  for j=0,nsite-1 do plot,rad[*,j],title='Actual rad. '+site_name[j],color=1
  device,/close
  set_plot,'x'
ENDIF

; load cloudiness
ncdf_read, cloud, file = cfile, /no_struct, variables = ['cloud']
cl = fltarr (nt, nsite)
FOR j = 0, nsite - 1 DO cl[*,j] = cloud[x[jsite[j]],y[jsite[j]],0:nt-1]
a = total(temporary(cloud))
IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  set_plot,'ps'
  device,file='cloud_sites.ps'
  for j=0,nsite-1 do plot,cl[*,j],title='Cloudiness '+site_name[j],color=1
  device,/close
  set_plot,'x'
ENDIF

; compute the statistics for the radiation ratio mean and range
FOR j = 0, nsite - 1 DO BEGIN

  ratio =  (rad[*,j] / rpot[*,j]) > 0 < 1.2

  FOR i=0,11 DO BEGIN

    l0 = where (p[*,j] GT 0 AND m EQ i)
    IF l0[0] EQ -1 THEN BEGIN
      mom1 = [missing, missing^2]
    ENDIF ELSE IF n_elements(l0) EQ 1 THEN BEGIN
      mom1 = [ratio[l0], missing^2]
    ENDIF ELSE BEGIN
      mom1 = moment (ratio[l0])
    ENDELSE

    ; daily mean cloudiness, mean and std. dev. on rain days
    mrmw[i,j] = mom1[0]
    mrsw[i,j] = sqrt(mom1[1])

    l0 = where (p[*,j] EQ 0 AND m EQ i)
    IF l0[0] EQ -1 THEN BEGIN
      mom1 = [missing, missing^2]
    ENDIF ELSE IF n_elements(l0) EQ 1 THEN BEGIN
      mom1 = [ratio[l0], missing^2]
    ENDIF ELSE BEGIN
      mom1 = moment (ratio[l0])
    ENDELSE

    ; daily mean cloudiness, mean and std. dev. on dry days
    mrmd[i,j] = mom1[0]
    mrsd[i,j] = sqrt(mom1[1])

    ; quadratic fit of ratio vs. cloudiness
    l0 = where (m EQ i)
    result = poly_fit (cl[l0,j], ratio[l0], 2)
    aratio[i,*,j] = result

    ; if missing value in rain or dry days, use from the other
    IF mrmw[i,j] EQ missing THEN mrmw[i,j] = mrmd[i,j]
    IF mrsw[i,j] EQ missing THEN mrsw[i,j] = mrsd[i,j]
    IF mrmd[i,j] EQ missing THEN mrmd[i,j] = mrmw[i,j]
    IF mrsd[i,j] EQ missing THEN mrsd[i,j] = mrsw[i,j]

  ENDFOR
  IF j mod 100 EQ 0 THEN print, 'Radiation, land point ', j, '   of ', nsite
ENDFOR

IF KEYWORD_SET (sites) THEN BEGIN
; plot the results
  set_plot,'ps'
  device,file='sites_weather_statistics.ps',/color
  FOR j = 0, nsite - 1 DO BEGIN
    plot,findgen(12)+1,mp[*,j],$
      title='Precip., wet prob. (thin), '+$
            ' wet-wet prob. (- - -), dry-wet prob. (...) '+ $
      site_name[j],thick=5,xrange=[1,12],/xstyle,xtitle='Month',$
      ytitle = 'mm/day OR fraction',color=1
    oplot,findgen(12)+1,mw[*,j],color=1
    oplot,findgen(12)+1,mww[*,j],li=2,color=1
    oplot,findgen(12)+1,mwd[*,j],li=1,color=1
  ENDFOR
  FOR j = 0, nsite - 1 DO BEGIN
    plot,findgen(12)+1,mtmd[*,j],$
      title='Temperature mean and SD (thin) (- - -: wet days) '+$
      site_name[j],thick=5,xrange=[1,12],/xstyle,xtitle='Month',$
      ytitle = 'deg C',color=1
    oplot,findgen(12)+1,mtsd[*,j],color=1
    oplot,findgen(12)+1,mtmw[*,j],li=2,thick=5,color=1
    oplot,findgen(12)+1,mtsw[*,j],li=2,color=1
  ENDFOR
  FOR j = 0, nsite - 1 DO BEGIN
    plot,findgen(12)+1,mtrmd[*,j],$
      title='Temperature range, mean and SD (thin) (- - -: wet days) '+$
      site_name[j],thick=5,xrange=[1,12],/xstyle,xtitle='Month',$
      ytitle = 'deg C',color=1
    oplot,findgen(12)+1,mtrsd[*,j],color=1
    oplot,findgen(12)+1,mtrmw[*,j],li=2,thick=5,color=1
    oplot,findgen(12)+1,mtrsw[*,j],li=2,color=1
  ENDFOR
  FOR j = 0, nsite - 1 DO BEGIN
    plot,findgen(12)+1,mrmd[*,j],color=1,$
      title='Rad. ratio mean and SD (thin) (- - -: wet days) '+$
      site_name[j],thick=5,xrange=[1,12],/xstyle,xtitle='Month'
    oplot,findgen(12)+1,mrsd[*,j],color=1
    oplot,findgen(12)+1,mrmw[*,j],li=2,thick=5,color=1
    oplot,findgen(12)+1,mrsw[*,j],li=2,color=1
  ENDFOR
  FOR j = 0, nsite - 1 DO BEGIN
    cld = findgen(101)/100
    plot,cld,aratio(0,0,j)+aratio(0,1,j)*cld+aratio(0,2,j)*cld^2,$
      yrange=[0,1],title='Rad. ratio to cloud fraction '+$
              ' Jan (black via blue green orange) to Dec  '+ $
        site_name[j],xtitle='Cloudiness [fraction]',$
        ytitle = 'actual to potential radiation', color=1
    for i=1,11 do oplot,cld,aratio(i,0,j)+aratio(i,1,j)*cld+$
      aratio(i,2,j)*cld^2,color=i*20
  ENDFOR
  device,/close
  set_plot,'x'
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (4) OUTPUT TO NETCDF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; output file name & graphics
tag = strcompress(year0_obs,/rem)+'-'+strcompress(year1_obs,/rem)+'.nc'
IF KEYWORD_SET (sites) THEN BEGIN
  wfile = dir + 'climate_pars_sites_'+tag; graphics
  IF KEYWORD_SET (lores) THEN BEGIN
    print, 'Combination of keywords SITES and LORES not allowed.'
    STOP
  ENDIF
ENDIF ELSE BEGIN
  wfile = dir + 'climate_pars_'+tag
  IF KEYWORD_SET (lores) THEN $
    wfile = dir + 'climate_pars_lores_'+tag
ENDELSE

print, 'Writing results to file...'

;save climate statistics in NetCDF file
id = NCDF_CREATE(wfile, /CLOBBER)

; define dimensions
dtime  = NCDF_DIMDEF(id, 'time', 12)
dorder  = NCDF_DIMDEF(id, 'order', 3)
dgrid   = NCDF_DIMDEF(id, 'gridnumber', /UNLIMITED)
IF KEYWORD_SET (sites) THEN BEGIN
  s = size (byte(site_name))
  dstrl = NCDF_DIMDEF(id, 'strl', s[1])
ENDIF

vtime = NCDF_VARDEF(id, 'time', [dtime], /SHORT)
NCDF_ATTPUT, id, vtime, 'long_name', 'month of year', /CHAR
NCDF_ATTPUT, id, vtime, 'units', 'month', /CHAR

; define variables
IF KEYWORD_SET (sites) THEN BEGIN
  vsname = NCDF_VARDEF(id, 'sitename', [dstrl,dgrid], /CHAR)
  NCDF_ATTPUT, id, vsname, 'long_name', 'name of site', /CHAR
ENDIF

vlon = NCDF_VARDEF(id, 'lon', [dgrid], /FLOAT)
NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

vlat = NCDF_VARDEF(id, 'lat', [dgrid], /FLOAT)
NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

vxcoord = NCDF_VARDEF(id, 'xcoord', [dgrid], /LONG)
NCDF_ATTPUT, id, vxcoord, 'long_name', 'integer horizontal grid coordinate', /CHAR

vycoord = NCDF_VARDEF(id, 'ycoord', [dgrid], /LONG)
NCDF_ATTPUT, id, vycoord, 'long_name', 'integer vertical grid coordinate', /CHAR

vmp = NCDF_VARDEF(id, 'mp', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmp, 'long_name', 'Mean monthly precipitation', /CHAR
NCDF_ATTPUT, id, vmp, 'units', 'mm/day', /CHAR
NCDF_ATTPUT, id, vmp, 'missing_value', smissing, /FLOAT
NCDF_ATTPUT, id, vmp, 'minimum_value', precip_min, /FLOAT


vmww = NCDF_VARDEF(id, 'mww', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmww, 'long_name', 'Mean monthly wet-to-wet day transition probability', /CHAR
NCDF_ATTPUT, id, vmww, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vmww, 'missing_value', smissing, /FLOAT

vmwd = NCDF_VARDEF(id, 'mwd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmwd, 'long_name', 'Mean monthly dry-to-wet day transition probability', /CHAR
NCDF_ATTPUT, id, vmwd, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vmwd, 'missing_value', smissing, /FLOAT

vmpkwb = NCDF_VARDEF(id, 'mpkwb', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmpkwb, 'long_name', 'Scale parameter of Weibull distribution for precipitation on rainy days', /CHAR
NCDF_ATTPUT, id, vmpkwb, 'units', 'unitless', /CHAR
NCDF_ATTPUT, id, vmpkwb, 'missing_value', smissing, /FLOAT

vmtmw = NCDF_VARDEF(id, 'mtmw', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtmw, 'long_name', 'Mean monthly temperature on rain days', /CHAR
NCDF_ATTPUT, id, vmtmw, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtmw, 'missing_value', smissing, /FLOAT

vmtsw = NCDF_VARDEF(id, 'mtsw', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtsw, 'long_name', 'Std. dev. of monthly temperature on rain days', /CHAR
NCDF_ATTPUT, id, vmtsw, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtsw, 'missing_value', smissing, /FLOAT

vmtrmw = NCDF_VARDEF(id, 'mtrmw', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtrmw, 'long_name', 'Mean monthly diurnal temperature range on rain days', /CHAR
NCDF_ATTPUT, id, vmtrmw, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtrmw, 'missing_value', smissing, /FLOAT

vmtrsw = NCDF_VARDEF(id, 'mtrsw', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtrsw, 'long_name', 'Std. dev. of monthly diurnal temperature range on rain days', /CHAR
NCDF_ATTPUT, id, vmtrsw, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtrsw, 'missing_value', smissing, /FLOAT

vmtmd = NCDF_VARDEF(id, 'mtmd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtmd, 'long_name', 'Mean monthly temperature on dry days', /CHAR
NCDF_ATTPUT, id, vmtmd, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtmd, 'missing_value', smissing, /FLOAT

vmtsd = NCDF_VARDEF(id, 'mtsd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtsd, 'long_name', 'Std. dev. of monthly temperature on dry days', /CHAR
NCDF_ATTPUT, id, vmtsd, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtsd, 'missing_value', smissing, /FLOAT

vmtrmd = NCDF_VARDEF(id, 'mtrmd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtrmd, 'long_name', 'Mean monthly diurnal temperature range on dry days', /CHAR
NCDF_ATTPUT, id, vmtrmd, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtrmd, 'missing_value', smissing, /FLOAT

vmtrsd = NCDF_VARDEF(id, 'mtrsd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmtrsd, 'long_name', 'Std. dev. of monthly diurnal temperature range on dry days', /CHAR
NCDF_ATTPUT, id, vmtrsd, 'units', 'degrees Celsius', /CHAR
NCDF_ATTPUT, id, vmtrsd, 'missing_value', smissing, /FLOAT

vmrmw = NCDF_VARDEF(id, 'mrmw', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmrmw, 'long_name', 'Mean monthly ratio of actual to potential solar incoming radiation on rain days', /CHAR
NCDF_ATTPUT, id, vmrmw, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vmrmw, 'missing_value', smissing, /FLOAT

vmrsw = NCDF_VARDEF(id, 'mrsw', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmrsw, 'long_name', 'Std. dev. of monthly ratio of actual to potential solar incoming radiation on rain days', /CHAR
NCDF_ATTPUT, id, vmrsw, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vmrsw, 'missing_value', smissing, /FLOAT

vmrmd = NCDF_VARDEF(id, 'mrmd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmrmd, 'long_name', 'Mean monthly ratio of actual to potential solar incoming radiation on dry days', /CHAR
NCDF_ATTPUT, id, vmrmd, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vmrmd, 'missing_value', smissing, /FLOAT

vmrsd = NCDF_VARDEF(id, 'mrsd', [dtime,dgrid], /FLOAT)
NCDF_ATTPUT, id, vmrsd, 'long_name', 'Std. dev. of monthly ratio of actual to potential solar incoming radiation on dry days', /CHAR
NCDF_ATTPUT, id, vmrsd, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vmrsd, 'missing_value', smissing, /FLOAT

varatio = NCDF_VARDEF(id, 'aratio', [dtime,dorder,dgrid], /FLOAT)
NCDF_ATTPUT, id, varatio, 'long_name', 'Parameters of quadratic fit of ratio of actual to potential solar incoming radiation to cloudiness', /CHAR
NCDF_ATTPUT, id, varatio, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, varatio, 'missing_value', smissing, /FLOAT

; define global attributes
NCDF_ATTPUT, id, /GLOBAL, 'title', 'climate statistics from VIC 1979 to 2005 data', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR

; change to DATA mode and put data into file
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, vtime, indgen(12) + 1
IF KEYWORD_SET (sites) THEN NCDF_VARPUT, id, vsname, byte(site_name)
NCDF_VARPUT, id, vlon, grid.lon[x[jsite]]
NCDF_VARPUT, id, vlat, grid.lat[y[jsite]]
IF NOT KEYWORD_SET (lores) THEN BEGIN
  NCDF_VARPUT, id, vxcoord, x[jsite]
  NCDF_VARPUT, id, vycoord, y[jsite]
ENDIF ELSE BEGIN
  NCDF_VARPUT, id, vxcoord, x0
  NCDF_VARPUT, id, vycoord, y0
ENDELSE
NCDF_VARPUT, id, vmp, mp
NCDF_VARPUT, id, vmww, mww
NCDF_VARPUT, id, vmwd, mwd
NCDF_VARPUT, id, vmpkwb, mpkwb
NCDF_VARPUT, id, vmtmw, mtmw
NCDF_VARPUT, id, vmtsw, mtsw
NCDF_VARPUT, id, vmtrmw, mtrmw
NCDF_VARPUT, id, vmtrsw, mtrsw
NCDF_VARPUT, id, vmtmd, mtmd
NCDF_VARPUT, id, vmtsd, mtsd
NCDF_VARPUT, id, vmtrmd, mtrmd
NCDF_VARPUT, id, vmtrsd, mtrsd
NCDF_VARPUT, id, vmrmw, mrmw
NCDF_VARPUT, id, vmrsw, mrsw
NCDF_VARPUT, id, vmrmd, mrmd
NCDF_VARPUT, id, vmrsd, mrsd
NCDF_VARPUT, id, varatio, aratio
; close file
NCDF_CLOSE, id

; checking the results
;device,decompose=0
;loadct,41
;
;ncdf_read,file='/Users/wolfgang/Models/CCDAS/input/climate/climate_pars_1979-2005.nc',data,/all,attr=attr
;help,/str,data
;n=n_elements(data.lon)
;d=fltarr(180,90,12)+1e6
;for j=0,n-1 do d[data.xcoord[j],data.ycoord[j],*]=data.mwd[*,j]
;for j=0,11 do begin & tvf,/bar,missing=1e6,d[*,*,j],title=j & wait,0.5 & endfor
;tvf,/bar,missing=1e6,d[*,*,6]
;tvf,/bar,missing=1e6,total(d,3)/12
;
;ncdf_read,file='/Users/wolfgang/Models/CCDAS/input/climate/climate_pars_lores_1979-2005.nc',data,/all,attr=attr
;d=fltarr(max(data.xcoord)+1,max(data.ycoord)+1,12)+1e6
;n=n_elements(data.lon)
;for j=0,n-1 do d[data.xcoord[j],data.ycoord[j],*]=data.mwd[*,j]
;for j=0,11 do begin & tvf,/bar,missing=1e6,d[*,*,j],title=j & wait,0.5 & endfor

END












