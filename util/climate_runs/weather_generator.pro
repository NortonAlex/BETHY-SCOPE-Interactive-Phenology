PRO weather_generator, sites = sites, lores = lores, $
                       gcm = gcm, scenario=scenario
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; WEATHER GENERATOR
;
; GENERATES DAILY WEATHER FROM year0 TO year1 FROM MONTHLY
; INPUT DATA
;
; UNTIL END OF OBSERVATIONAL RECORD, THE MONTHLY DATA USED ARE FROM OBSERVATIONS
; AFTER THAT YEAR, ANOMALY CORRECTED CLIMATE SCENARIOS
;
; NEEDS SEVERAL PREVIOUS CONVERSION STEPS TO EXTRACT MEANS AND PARAMETERS
; SEE ACCOMPANYING BATCH FILE scenario.bat
;
; AUTHOR:
; Wolfgang Knorr, March-April 2008
; Update in August 2009
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (1) Declarations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
;dir = '/Users/wolfgang/Models/CCDAS/input/climate/'
;sdir = '/Users/wolfgang/Models/CCDAS/input/eddy_sites/'
dir = '/Volumes/Work/ccdas_input/climate/'
sdir = '/Volumes/Work/ccdas_input/eddy_sites/'
missing = 1e6
smissing = '1e6'
year0 = 1979  ; start year for weather generator
year1 = 2039  ; end year for weather generator
precip_min = 0.01 ; minimum precipitation on wet day [mm]
precip_max = 1000. ; maximum daily precipitation [mm]
precip_mmax = 10000. ; maximum monthly precipitation [mm]
temp_max = 60.0 ; maximum temperature [degC]
temp_min = -80.0 ; minimum temperature [degC]
IF NOT KEYWORD_SET (gcm) THEN gcm = 'HadCM3'
IF NOT KEYWORD_SET (scenario) THEN scenario = 'A1'
out_range = strcompress(year0,/rem)+'_'+strcompress(year1,/rem)
restag = 'hires'
IF KEYWORD_SET (lores) THEN restag = 'lores'

; input files
file_scen = dir + gcm + '_' + scenario + '_2008_2099_monthly.nc'
file_clim = dir + gcm + '_' + scenario+'_1979_2007_monthly_climatology.nc'
file_obs = dir + 'forcing.hires.monthly.1979-2008.nc'
IF KEYWORD_SET (sites) THEN BEGIN
  file_para = dir + 'climate_pars_sites_1979-2007.nc'
ENDIF ELSE IF KEYWORD_SET (lores) THEN BEGIN
  file_para = dir + 'climate_pars_lores_1979-2007.nc'
ENDIF ELSE BEGIN
  file_para = dir + 'climate_pars_1979-2007.nc'
ENDELSE

; output file names for default option (global)
wfile_precip = dir+gcm+'_'+scenario+'_precip_daily_'+restag+'_'+out_range+'.nc'
wfile_tmin = dir+gcm+'_'+scenario+'_tmin_daily_'+restag+'_'+out_range+'.nc'
wfile_tmax = dir+gcm+'_'+scenario+'_tmax_daily_'+restag+'_'+out_range+'.nc'
wfile_swdown = dir+gcm+'_'+scenario+'_swdown_daily_'+restag+'_'+out_range+'.nc'

; output file name element for option 'sites'
sfile = sdir + gcm + '_' + scenario + '_daily_'+out_range+'_'

; calendar stuff
dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (2) Initialization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ncdf_read, file = file_scen, scen, /all
ncdf_read, file = file_clim, clim, /all
ncdf_read, file = file_para, para, /all
ncdf_read, file = file_obs, obs, /all

; load elevation
ncdf_read, file = $
;    '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', $
    '/Volumes/Work/wolfgang/CCDAS/ppm/control_bethy/bethy_grid3veg_v2.nc', $
    /no_struct, variables = ['elev'], elev

; grid information (or list of site locations)
lon = para.lon
lat = para.lat
; x,y coordinates on the specified grid
x0 = para.xcoord
y0 = para.ycoord
IF KEYWORD_SET (lores) THEN BEGIN
  ; hard-wired 2x2 degree grid
  nx = 180
  ny = 90
  ; x,y coordinates mapped to 2x2 degrees
  x = fix((lon+180.)/360.*nx)
  y = fix((lat+90.)/180.*ny)
ENDIF ELSE BEGIN
  x = x0
  y = y0
ENDELSE
nsite = n_elements(lon)
IF KEYWORD_SET (sites) THEN site_name = strcompress(para.sitename,/rem)

; compute pointers to months  
l1 = where (obs.time ge year0)
m0 = l1[0] ; pointer to first month in observations
IF fix(obs.time[m0]) NE year0 THEN STOP
m1 = m0 + n_elements(obs.time) - 1 ; pointer to last month in observations
l1 = where (scen.time gt obs.time[m1]+0.5/12)
m2 = l1[0] ; pointer to first month in scenario
l1 = where (scen.time lt year1+1)
m3 = l1[n_elements(l1)-1]
IF fix(scen.time[m3]) NE year1 THEN STOP
print, 'Observational period (monthly): ',obs.time[m0],obs.time[m1]
print, 'Scenario period (monthly): ',scen.time[m2],scen.time[m3]
; check if gap exists between observations and scenario
IF scen.time[m2] - obs.time[m1] GT 1.001/12. THEN STOP
nyr = year1 - year0 + 1
nm = nyr * 12L
time = fltarr (nm)
time[0:m1-m0] = obs.time[m0:m1]
time[m1-m0+1:nm-1] = scen.time[m2:m3]
lm0 = fix(obs.time[m0:m1]*12 mod 12) ; pointer to month in year for obs. period
lm2 = fix(scen.time[m2:m3]*12 mod 12) ; pointer to month in year for scen. period
nyr0 = intarr(12)
for i = 0, m1-m0 do nyr0[lm0[i]] = nyr0[lm0[i]] + 1 ; number of years per month in obs. period

; compute length of time series in days
n_leap = fix(total((fix(time) mod 4) eq 0)/12)
nt = nyr*365L + n_leap

; generate calendar data (day, month, year)
yr = intarr(nt) & month = intarr(nt) & m = intarr(nt) & d = intarr(nt) & jday = intarr(nt)
mth = intarr(nm)
d[0] = 0 & month[0] = 0 & m[0] = 0 & yr[0] = 0 & jday[0] = 1 & mth[0] = 0
FOR i=1,nt-1 DO BEGIN
  d[i] = d[i-1] + 1
  jday[i] = jday[i-1] + 1
  month[i] = month[i-1]
  m[i] = m[i-1]
  yr[i] = yr[i-1]
  dpm0 = dpm[month[i]]
  IF month[i] EQ 1 AND ((yr[i]+year0) MOD 4) EQ 0 THEN dpm0 = dpm0 + 1
  IF d[i] EQ dpm0 THEN BEGIN
    d[i] = 0
    IF month[i] EQ 11 THEN BEGIN
      yr[i] = yr[i] + 1
      jday[i] = 1
    ENDIF
    month[i] = (month[i] + 1) MOD 12
    m[i] = m[i] + 1
    mth[m[i]] = month[i]
  ENDIF
ENDFOR
if max(yr)+1 NE nyr THEN STOP

; the monthly fields
mp = fltarr (nm, nsite)  ; monthly precipitation [mm/day]
mt = fltarr (nm, nsite) ; monthly temperature [deg C]
mtr = fltarr (nm, nsite) ; monthly diurnal temperature range [deg C]
mr = fltarr (nm, nsite) ; monthly ratio of actual over potential SW rad.

; the monthly observed climatology fields
mp_obs_clim = fltarr (12, nsite)
mt_obs_clim = fltarr (12, nsite)
mtr_obs_clim = fltarr (12, nsite)
mr_obs_clim = fltarr (12, nsite)

; the daily (output) arrays
p = fltarr (nt, nsite)
tmean = fltarr (nt, nsite)
trange = fltarr (nt, nsite)
ratio = fltarr (nt, nsite)
rpot = fltarr (nt, nsite)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (3) construct monthly mean scenario (before daily)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; process quadratic functions of actual to potential radiation to cloud fraction
; exclude rad. ratio to cloud function with positive gradient at either end
; and replace by average function (except for constant element)
; FROM GLOBAL RUN:
; Average aratio:      0.854844   -0.0677555    -0.497829
aratio_av = [0.854844, -0.0677555, -0.497829]
para.aratio[*,0,*] = para.aratio[*,0,*] > 0. ; exclude negative values
mask = (para.aratio[*,1,*] LE 0.) AND $
       (para.aratio[*,1,*]+2.*para.aratio[*,2,*] LE 0.)
para.aratio[*,1,*] = para.aratio[*,1,*]*mask + aratio_av[1]*(mask EQ 0)
para.aratio[*,2,*] = para.aratio[*,2,*]*mask + aratio_av[2]*(mask EQ 0)

; monthly mean precipitation
for j = 0, nsite-1 do mp[0:m1-m0,j] = obs.precip[x[j],y[j],m0:m1]
for i = 0, m1-m0 do mp_obs_clim[lm0[i],*] = mp_obs_clim[lm0[i],*] + mp[m0+i,*] / nyr0[lm0[i]]
for j = 0, nsite-1 do mp[m1-m0+1:nm-1,j] = scen.precip[x[j],y[j],m2:m3] $
                                         / clim.precip[x[j],y[j],lm2] $
                                         * mp_obs_clim[lm2,j]
mp = mp < precip_mmax

; monthly mean temperature
for j = 0, nsite-1 do mt[0:m1-m0,j] = (obs.tmax[x[j],y[j],m0:m1] + $
                                       obs.tmin[x[j],y[j],m0:m1]) / 2.
for i = 0, m1-m0 do mt_obs_clim[lm0[i],*] = mt_obs_clim[lm0[i],*] + mt[m0+i,*] / nyr0[lm0[i]]
for j = 0, nsite-1 do mt[m1-m0+1:nm-1,j] = scen.temp[x[j],y[j],m2:m3] $
                                         - clim.temp[x[j],y[j],lm2] $
                                         + mt_obs_clim[lm2,j]

; monthly mean diurnal temperature range
for j = 0, nsite-1 do mtr[0:m1-m0,j] = (obs.tmax[x[j],y[j],m0:m1] - $
                                        obs.tmin[x[j],y[j],m0:m1])
for i = 0, m1-m0 do mtr_obs_clim[lm0[i],*] = mtr_obs_clim[lm0[i],*] + mtr[m0+i,*] / nyr0[lm0[i]]
for j = 0, nsite-1 do mtr[m1-m0+1:nm-1,j] = mtr_obs_clim[lm2,j]

; monthly mean actual over potential shortwave downwelling radiation
for j = 0, nsite-1 do mr[0:m1-m0,j] = para.aratio[lm0,0,j] + $
  para.aratio[lm0,1,j] * obs.cloud[x[j],y[j],m0:m1] + $
  para.aratio[lm0,2,j] * obs.cloud[x[j],y[j],m0:m1]^2 > 0.
for i = 0, m1-m0 do mr_obs_clim[lm0[i],*] = mr_obs_clim[lm0[i],*] + mr[m0+i,*] / nyr0[lm0[i]]
for j = 0, nsite-1 do mr[m1-m0+1:nm-1,j] = para.aratio[lm2,1,j] * $
     (scen.cloud[x[j],y[j],m2:m3] - clim.cloud[x[j],y[j],lm2]) + $
     para.aratio[lm2,2,j] * $
     (scen.cloud[x[j],y[j],m2:m3]^2 - clim.cloud[x[j],y[j],lm2]^2) + $
     mr_obs_clim[lm2,j] > 0.

IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  set_plot,'ps'
  device,file=gcm+'_'+scenario+'_monthly_sites.ps'
  FOR j = 0, nsite-1 DO BEGIN
    plot,time,/xst,mp[*,j],ytitle='mm/day',$
      title=site_name[j]+' monthly mean precipitation + linear trend'
    linreg,time,mp[*,j],a,b
    oplot,[year0,year1],[a+b*year0,a+b*year1]
    plot,time,/xst,mt[*,j],ytitle='deg C',$
      title=site_name[j]+' monthly mean temperature + linear trend'
    linreg,time,mt[*,j],a,b
    oplot,[year0,year1],[a+b*year0,a+b*year1]
    plot,time,/xst,mtr[*,j],ytitle='deg C',$
      title=site_name[j]+' diurnal temperature range'
    plot,time,/xst,mr[*,j],title=site_name[j]+$
      ' monthly actual / pot. SW radiation + linear trend'
    linreg,time,mr[*,j],a,b
    oplot,[year0,year1],[a+b*year0,a+b*year1]
  ENDFOR
  device,/close
  set_plot,'x'
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (4) generate daily weather
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters for radiation model
laps = 0.006
p0 = 1.01325e5
a0 = [1.00011, 0.34221E-1, 0.128E-2, 0.719E-3, 0.77E-4]
s0 = 1360.    ;  total radiation
b0 = 0.185
d0 = 0.4

ndays0 = 365.

; compute mean monthly wet-day fraction
mw = para.mwd/(1.-para.mww+para.mwd)

; generate daily rainfall
FOR i = 0, nt - 1 DO BEGIN

 IF i EQ 0 THEN BEGIN
    ; first day: set rain/dry according to mean wet-day occurrence
    p[i,*] = randomu (seed, nsite) LE mw[month[i],*]
  ENDIF ELSE BEGIN
    ; all other days: set rain/dry from transitional probabilities
    lw = where (p[i-1,*] GT 0.) & nw = n_elements(lw)
    ld = where (p[i-1,*] EQ 0.) & nd = n_elements(ld)
    if lw[0] ne -1 then p[i,lw] = randomu (seed, nw) LE para.mww[month[i],lw]
    if ld[0] ne -1 then p[i,ld] = randomu (seed, nd) LE para.mwd[month[i],ld]
  ENDELSE

  ; generate rain amounts from Weibull distribution
  lw = where (p[i,*]) & nw = n_elements(lw)
  IF lw[0] NE -1 THEN BEGIN
    kwb = para.mpkwb[month[i],*]
    l0 = where (kwb EQ missing)
    IF l0[0] NE -1 THEN kwb[l0] = 0.8 ; default Weibull parameter
    lwb = (mp[m[i],lw]-precip_min)/gamma(1.+1./kwb)
    ;for count=0,999 do $
    ; print,transpose(lwb*(-alog(randomu (seed, nw)))^(1./kwb) + precip_min)
    p[i,lw] = (lwb*(-alog(randomu (seed, nw)))^(1./kwb) > 0.) + precip_min
  ENDIF

  IF i mod 1000 EQ 0 THEN print, 'Precipitation day ', i, '    of ', nt

ENDFOR

; generate daily temperature
FOR i = 0, nt - 1 DO BEGIN

  ; generate diurnal mean and range (rain and dry days)
  lw = where (p[i,*]) & nw = n_elements(lw)
  ld = where (p[i,*] EQ 0.) & nd = n_elements(ld)
  IF lw[0] NE -1 THEN BEGIN
    tmean[i,lw] = para.mtmw[month[i],lw] $
                  + mt[m[i],lw] - mt_obs_clim[month[i],lw] $
                  + randomn (seed, nw) * para.mtsw[month[i],lw]
    trange[i,lw] = para.mtrmw[month[i],lw] $
                  + mtr[m[i],lw] - mtr_obs_clim[month[i],lw] $
                  + randomn (seed, nw) * para.mtrsw[month[i],lw]
  ENDIF
  IF ld[0] NE -1 THEN BEGIN
    tmean[i,ld] = para.mtmw[month[i],ld] $
                  + mt[m[i],ld] - mt_obs_clim[month[i],ld] $
                  + randomn (seed, nd) * para.mtsw[month[i],ld]
    trange[i,ld] = para.mtrmw[month[i],ld] $
                  + mtr[m[i],ld] - mtr_obs_clim[month[i],ld] $
                  + randomn (seed, nd) * para.mtrsw[month[i],ld]
  ENDIF

  IF i mod 1000 EQ 0 THEN print, 'Temperature day ', i, '    of ', nt

ENDFOR

; generate daily shortwave downwelling radiation
FOR i = 0, nt - 1 DO BEGIN

  ; compute potential radiation
  alpha1 = 2. * !pi * (jday[i] - 1) / ndays0
  alpha2 = alpha1 * 2.0

  delta= -23.4 * COS (2. * !pi * (jday[i] + 10) / ndays0 )

  spd = SIN (lat * !pi / 180.) * SIN (delta * !pi / 180.)
  cpd = COS (lat * !pi / 180.) * COS (delta * !pi / 180.)
  dbodsq = a0(0) + a0(1) * COS (alpha1) + a0(2) * SIN (alpha1) + $
           a0(3) * COS (alpha2) + a0(4) * SIN (alpha2)
  pair = p0*(1./(1.+elev[x,y]*laps/(tmean[i,*]/2.+273.15)))
  FOR hour = 0, 23 DO BEGIN
    h = hour * !pi / 12.
    coszen =   (spd - cpd * COS (h)) > 1e-12
    rtop = s0 * dbodsq * coszen
    tdir = EXP (-b0 / coszen * pair / p0)

    rpot[i,*] = (rpot[i,*] + rtop * (d0 + (1. - d0) * tdir) / 24.) > 0.
  ENDFOR

  ; generate radiation data
  lw = where (p[i,*]) & nw = n_elements(lw)
  ld = where (p[i,*] EQ 0.) & nd = n_elements(ld)
  if lw[0] NE -1 then $
    ratio[i,lw] = para.mrmw[month[i],lw] $
                  + mr[m[i],lw] - mr_obs_clim[month[i],lw] $
                  + randomn (seed, nw) * para.mrsw[month[i],lw]
  if ld[0] NE -1 then $
    ratio[i,ld] = para.mrmd[month[i],ld]  $
                  + mr[m[i],ld] - mr_obs_clim[month[i],ld] $
                  + randomn (seed, nd) * para.mrsw[month[i],ld]

  IF i mod 1000 EQ 0 THEN print, 'Radiation day ', i, '    of ', nt

ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (8) Readjustment to match monthly means
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; re-scale precipitation to prescribed monthly means
mp_gen = fltarr (nm, nsite)
nday = intarr (nm)
for i = 0, nm-1 do nday[i] = total(m EQ i)
for i = 0, nm-1 do mp_gen[i,*] = total(p[where(m EQ i),*],1) / nday[i]
for i = 0, nt-1 do p[i,*] = p[i,*] / (mp_gen[m[i],*] $
  + (mp_gen[m[i],*] EQ 0.)) * mp[m[i],*]

; adjust temperature diurnal mean and range to prescribed monthly means
mt_gen = fltarr (nm, nsite)
mtr_gen = fltarr (nm, nsite)
for i = 0, nm-1 do mt_gen[i,*] = total(tmean[where(m EQ i),*],1) / nday[i]
for i = 0, nm-1 do mtr_gen[i,*] = total(trange[where(m EQ i),*],1) / nday[i]
for i = 0, nt-1 do tmean[i,*] = tmean[i,*] - mt_gen[m[i],*] + mt[m[i],*]
for i = 0, nt-1 do trange[i,*] = trange[i,*] - mtr_gen[m[i],*] + mtr[m[i],*]
trange = trange > 0

; adjust pot. to actual SW radation ratio to prescribed monthly means
;ratio = ratio > 0 < 1
mr_gen = fltarr (nm, nsite)
for i = 0, nm-1 do mr_gen[i,*] = total(ratio[where(m EQ i),*],1) / nday[i]
for i = 0, nt-1 do ratio[i,*] = ratio[i,*] - mr_gen[m[i],*] + mr[m[i],*]
ratio = ratio > 0 < 1

IF KEYWORD_SET (sites) THEN BEGIN
  ; checking:
  fname = '_'+gcm+'_'+scenario+'_'+strcompress(year0,/rem)+'_'+strcompress(year1,/rem)+'.ps'
  set_plot,'ps'
  device,file='precip'+fname
  for j=0,nsite-1 do plot,(jday-1)/365.+yr+year0,p[*,j],$
    title=site_name[j]+' daily precipitation',/xst,ytitle='mm/day'
  mp0=fltarr(nm,nsite)
  for i=0,nm-1 do for j=0,nsite-1 do mp0[i,j]=mean(p[where(m eq i),j])
  for j=0,nsite-1 do begin & plot,time,mp0[*,j],title=site_name[j]+$
    ' monthly precipitation',/xst,ytitle='mm/day' $
    & oplot,time,mp[*,j],li=2 & endfor
  mw0=fltarr(12,nsite)
  for i=0,11 do for j=0,nsite-1 do mw0[i,j]=mean(p[where(month eq i),j] gt 0.)
  for j=0,nsite-1 do begin & plot,indgen(12)+1,mw0[*,j],$
    /xst,xtitle='month',title=site_name[j]+' wet day fraction' & $
    oplot,indgen(12)+1,mw[*,j],li=2 & endfor
  device,/close
  device,file='temp'+fname
  for j=0,nsite-1 do plot,(jday-1)/365.+yr+year0,tmean[*,j],$
    title=site_name[j]+' daily mean temperature',/xst,ytitle='deg C'
  for j=0,nsite-1 do plot,(jday-1)/365.+yr+year0,trange[*,j],$
    title=site_name[j]+' daily temperature range',/xst,ytitle='deg C'
  mt0=fltarr(nm,nsite)
  for i=0,nm-1 do for j=0,nsite-1 do mt0[i,j]=mean(tmean[where(m eq i),j])
  for j=0,nsite-1 do begin & plot,time,mt0[*,j],title=site_name[j]+$
    ' monthly mean temperature',/xst,ytitle='deg C' $
    & oplot,time,mt[*,j],li=2 & endfor
  device,/close
  device,file='swdown'+fname
  for j=0,nsite-1 do plot,(jday-1)/365.+yr+year0,ratio[*,j]*rpot[*,j],$
    title=site_name[j]+' daily mean SW downw. rad.',/xst,ytitle='W/m^2'
  mswdown0=fltarr(nm,nsite)
  for i=0,nm-1 do for j=0,nsite-1 do begin & l1=where(m eq i) & $
    mswdown0[i,j]=mean(ratio[l1,j]*rpot[l1,j]) & endfor
  mrpot0=fltarr(nm,nsite)
  for i=0,nm-1 do for j=0,nsite-1 do mrpot0[i,j]=mean(rpot[where(m eq i),j])
  for j=0,nsite-1 do begin & plot,time,mswdown0[*,j],title=site_name[j]+$
    ' monthly mean SW downw. rad.',/xst,ytitle='W/m^2' $
    & oplot,time,mrpot0[*,j]*mr[*,j],li=2 & endfor
  device,/close
  set_plot,'x'
ENDIF

print, 'Saving to files...'

IF KEYWORD_SET (sites) THEN BEGIN

  ; save in one file per site
  FOR j = 0, nsite-1 DO BEGIN

    id = NCDF_CREATE(sfile+site_name[j]+'.nc', /CLOBBER)

    ; define dimensions
    dtime  = NCDF_DIMDEF(id, 'time', nt)

    vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
    NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' + $
      strcompress(year0,/rem), /CHAR
    NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR

    ; define variables
    vlon = NCDF_VARDEF(id, 'lon', /FLOAT)
    NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
    NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

    vlat = NCDF_VARDEF(id, 'lat', /FLOAT)
    NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
    NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

    vxcoord = NCDF_VARDEF(id, 'xcoord', /LONG)
    NCDF_ATTPUT, id, vxcoord, 'long_name', $
      'integer horizontal grid coordinate', /CHAR

    vycoord = NCDF_VARDEF(id, 'ycoord', /LONG)
    NCDF_ATTPUT, id, vycoord, 'long_name', $
      'integer vertical grid coordinate', /CHAR

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
      + strcompress(year0,/rem) + ' to ' + strcompress(year1,/rem), /CHAR
    NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
    NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
    NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR

    ; change to DATA mode and put data into file
    NCDF_CONTROL, id, /ENDEF

    NCDF_VARPUT, id, vlon, lon[j]
    NCDF_VARPUT, id, vlat, lat[j]
    NCDF_VARPUT, id, vxcoord, x[j]
    NCDF_VARPUT, id, vycoord, y[j]
    NCDF_VARPUT, id, vtime, indgen(nt)
    NCDF_VARPUT, id, vprecip, p[*,j] > 0. < precip_max
    NCDF_VARPUT, id, vtmin, (tmean[*,j] - trange[*,j] / 2.) $
                     < temp_max > temp_min
    NCDF_VARPUT, id, vtmax, (tmean[*,j] + trange[*,j] / 2.) $
                     < temp_max > temp_min
    NCDF_VARPUT, id, vswdown, ratio[*,j] * rpot[*,j]
    ; close file
    NCDF_CLOSE, id

  ENDFOR

ENDIF ELSE BEGIN

  IF KEYWORD_SET (lores) THEN BEGIN
    ; save not the co-ordinates on 2x2 grid, but on lores grid
    x = x0
    y = y0
  ENDIF

  ;save precipitation data in NetCDF files
  id = NCDF_CREATE(wfile_precip, /CLOBBER)
  
  ; define dimensions
  dgrid   = NCDF_DIMDEF(id, 'gridnumber', nsite)
  dtime  = NCDF_DIMDEF(id, 'time', nt)
  
  vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
  NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' $
    + strcompress(year0,/rem), /CHAR
  NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR
  
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
  
  vprecip = NCDF_VARDEF(id, 'precip', [dtime, dgrid], /FLOAT)
  NCDF_ATTPUT, id, vprecip, 'long_name', 'Daily precipitation', /CHAR
  NCDF_ATTPUT, id, vprecip, 'units', 'mm', /CHAR
  NCDF_ATTPUT, id, vprecip, 'missing_value', smissing, /FLOAT
  
  ; define global attributes
  NCDF_ATTPUT, id, /GLOBAL, 'title', $
    'daily generated climate data for year ' + $
    strcompress(year0,/rem) + ' to ' + strcompress(year1,/rem), /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR
  
  ; change to DATA mode and put data into file
  NCDF_CONTROL, id, /ENDEF
  
  NCDF_VARPUT, id, vlon, lon
  NCDF_VARPUT, id, vlat, lat
  NCDF_VARPUT, id, vxcoord, x
  NCDF_VARPUT, id, vycoord, y
  NCDF_VARPUT, id, vtime, indgen(nt)
  NCDF_VARPUT, id, vprecip, p > 0. < precip_max
  ; close file
  NCDF_CLOSE, id

  print, 'Number of cases with P < 0: ', $
    fix(total(p LT 0))
  print, 'Number of cases with P > ', precip_max, ': ', $
    fix(total(p GT precip_max))

  ;save minimum temperature data in NetCDF files
  id = NCDF_CREATE(wfile_tmin, /CLOBBER)
  
  ; define dimensions
  dgrid   = NCDF_DIMDEF(id, 'gridnumber', nsite)
  dtime  = NCDF_DIMDEF(id, 'time', nt)
  
  vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
  NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' + $
    strcompress(year0,/rem), /CHAR
  NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR
  
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
  
  vtmin = NCDF_VARDEF(id, 'tmin', [dtime, dgrid], /FLOAT)
  NCDF_ATTPUT, id, vtmin, 'long_name', 'Daily minimum temperatue', /CHAR
  NCDF_ATTPUT, id, vtmin, 'units', 'degrees Celsius', /CHAR
  NCDF_ATTPUT, id, vtmin, 'missing_value', smissing, /FLOAT
  
  ; define global attributes
  NCDF_ATTPUT, id, /GLOBAL, 'title', $
    'daily generated climate data for year ' + $
    strcompress(year0,/rem) + ' to ' + strcompress(year1,/rem), /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR
  
  ; change to DATA mode and put data into file
  NCDF_CONTROL, id, /ENDEF
  
  NCDF_VARPUT, id, vlon, lon
  NCDF_VARPUT, id, vlat, lat
  NCDF_VARPUT, id, vxcoord, x
  NCDF_VARPUT, id, vycoord, y
  NCDF_VARPUT, id, vtime, indgen(nt)
  NCDF_VARPUT, id, vtmin, (tmean-trange/2.) < temp_max > temp_min
  ; close file
  NCDF_CLOSE, id

  print, 'Number of cases with Tmin < ', temp_min, ': ', $
    fix(total((tmean-trange/2.) LT temp_min))
  print, 'Number of cases with Tmin > ', temp_max, ': ', $
    fix(total((tmean-trange/2.) GT temp_max))
  
  ;save climate statistics in NetCDF files
  id = NCDF_CREATE(wfile_tmax, /CLOBBER)
  
  ; define dimensions
  dgrid   = NCDF_DIMDEF(id, 'gridnumber', nsite)
  dtime  = NCDF_DIMDEF(id, 'time', nt)
  
  vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
  NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' $
    + strcompress(year0,/rem), /CHAR
  NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR
  
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
  
  vtmax = NCDF_VARDEF(id, 'tmax', [dtime, dgrid], /FLOAT)
  NCDF_ATTPUT, id, vtmax, 'long_name', 'Daily maximum temperatue', /CHAR
  NCDF_ATTPUT, id, vtmax, 'units', 'degrees Celsius', /CHAR
  NCDF_ATTPUT, id, vtmax, 'missing_value', smissing, /FLOAT
  
  ; define global attributes
  NCDF_ATTPUT, id, /GLOBAL, 'title', $
    'daily generated climate data for year ' + $
    strcompress(year0,/rem) + ' to ' + strcompress(year1,/rem), /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR
  
  ; change to DATA mode and put data into file
  NCDF_CONTROL, id, /ENDEF
  
  NCDF_VARPUT, id, vlon, lon
  NCDF_VARPUT, id, vlat, lat
  NCDF_VARPUT, id, vxcoord, x
  NCDF_VARPUT, id, vycoord, y
  NCDF_VARPUT, id, vtime, indgen(nt)
  NCDF_VARPUT, id, vtmax, (tmean+trange/2.) < temp_max > temp_min
  ; close file
  NCDF_CLOSE, id
  
  print, 'Number of cases with Tmax < ', temp_min, ': ', $
    fix(total((tmean+trange/2.) LT temp_min))
  print, 'Number of cases with Tmax > ', temp_max, ': ', $
    fix(total((tmean+trange/2.) GT temp_max))
  
  ;save climate statistics in NetCDF file
  id = NCDF_CREATE(wfile_swdown, /CLOBBER)
  
  ; define dimensions
  dgrid   = NCDF_DIMDEF(id, 'gridnumber', nsite)
  dtime  = NCDF_DIMDEF(id, 'time', nt)
  
  vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
  NCDF_ATTPUT, id, vtime, 'long_name', 'days since 1 Jan ' + $
    strcompress(year0,/rem), /CHAR
  NCDF_ATTPUT, id, vtime, 'units', 'days', /CHAR
  
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
  
  vswdown = NCDF_VARDEF(id, 'swdown', [dtime, dgrid], /FLOAT)
  NCDF_ATTPUT, id, vswdown, 'long_name', $
    'Daily mean incoming solar radation', /CHAR
  NCDF_ATTPUT, id, vswdown, 'units', 'W / m^2', /CHAR
  NCDF_ATTPUT, id, vswdown, 'missing_value', smissing, /FLOAT
  
  ; define global attributes
  NCDF_ATTPUT, id, /GLOBAL, 'title', $
    'daily generated climate data for year ' + $
    strcompress(year0,/rem) + ' to ' + strcompress(year1,/rem), /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
  NCDF_ATTPUT, id, /GLOBAL, 'creation_date', 'March 2008', /CHAR
  
  ; change to DATA mode and put data into file
  NCDF_CONTROL, id, /ENDEF
  
  NCDF_VARPUT, id, vlon, lon
  NCDF_VARPUT, id, vlat, lat
  NCDF_VARPUT, id, vxcoord, x
  NCDF_VARPUT, id, vycoord, y
  NCDF_VARPUT, id, vtime, indgen(nt)
  NCDF_VARPUT, id, vswdown, ratio * rpot
  ; close file
  NCDF_CLOSE, id

ENDELSE

; checking the results
;device,decompose=0
;loadct,41
;ncdf_read,file='/Users/wolfgang/Models/CCDAS/input/climate/ECHAM5_A1_precip_daily_hires_1979_2039.nc',data,/all,attr=attr
;ncdf_read,file='/Users/wolfgang/Models/CCDAS/input/climate/ECHAM5_A1_precip_daily_lores_1979_2039.nc',data,/all,attr=attr
;help,/str,data
;n=n_elements(data.lon)
;nt=n_elements(data.time)
;j=n/2
;print,data.lon[j],data.lat[j]
;plot,data.precip[*,j]
;d=fltarr(max(data.xcoord)+1,max(data.ycoord)+1)+1e6
;FOR t=nt*0.75,nt-1 DO BEGIN
;  for j=0,n-1 do d[data.xcoord[j],data.ycoord[j]]=data.precip[t,j]
;  tvf,/bar,missing=1e6,d,min=0,max=100.,title='Day '+string(t)
;  wait,0.2
;ENDFOR

END

