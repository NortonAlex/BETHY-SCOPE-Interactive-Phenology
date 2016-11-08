PRO fapar_compare, ps = ps, site_file = site_file, prior = prior, pft_plot = pft_plot

; PURPOSE:
; plots simulated and observed FAPAR, observed with error bars
; OPTIONS:
; ps         generates postcript in file 'fapar_sites.ps' in current dir.
;            for ps = 2, error bars are plotted in grey
; site_file  name of site specification file, which is looked for
;            in dir. ../control_bethy (so running this out of /util
;            will work fine).
; prior      additionally looks for results in the directory /output_prior
;            instead of /output for results of prior (uncalibrated) runs
;            created by running 'post' before optimisation and
;            copying 'cp output output_prior'
; pft_plot   also plot contribution of invidual PFTs

; read sites information
IF NOT KEYWORD_SET(site_file) THEN $
  site_file = 'site_specs_all.dat'
ng = 0
nsite = 0
nsmax = 10000
maxpft = 3
pft = intarr(maxpft,nsmax)
frac = fltarr(maxpft,nsmax)
site = strarr(nsmax)
coord = intarr(4,nsmax)
pftt = intarr(maxpft)
fract = fltarr(maxpft)
coordt = intarr(4)
openr,u,/get_lun,'../control_bethy/'+site_file
line = ' '
readf,u,line
WHILE NOT EOF(u) DO BEGIN
  readf,u,line
  reads,line,lat,lon,elev,pftt,fract,coordt
  print,pftt,fract,coordt
  n = fix(total(pftt NE 0)) ; ignore sub-gridcells with PFT=0
  pft[*,nsite] = pftt
  frac[*,nsite] = fract
  coord[*,nsite] = coordt
  site[nsite] = strcompress(strmid(line,strpos(line,'!')+1),/rem)
  ng = ng + n
  nsite = nsite + 1
ENDWHILE
free_lun, u
pft = pft[*,0:nsite-1] ; the PFT code
frac = frac[*,0:nsite-1] ; the fraction of this PFT within the sub-site
coord = coord[*,0:nsite-1] ; the coordinates within the scene
site = site[0:nsite-1] ; the site name (may be the same for several sub-sites)

; MODELLED FAPAR
dir='../input/eddy_sites/'
mdir='../output/'
mdir0='../output_prior/'
openr,u,mdir+'fpar1_site.dat',/get_lun
readf,u,lat,lon,nyr,year0
free_lun,u
mfpar=fltarr(365,nyr,maxpft,nsite)
mfpar0=fltarr(365,nyr,maxpft,nsite)
yr0=2
yr1=3
in=fltarr(365)
;
FOR k=0,2 DO BEGIN
  openr,u,mdir+'fpar'+string(k+1,fo='(I1)')+'_site.dat',/get_lun
  FOR s=0,nsite-1 DO BEGIN
    readf,u,lat,lon,nyr,year0
    FOR i=0,nyr-1 DO BEGIN
      readf,u,in & mfpar[*,i,k,s]=in
    ENDFOR
  ENDFOR
  free_lun,u
ENDFOR
mfpar=reform(mfpar[*,yr0:yr1,*,*],(yr1-yr0+1)*365,maxpft,nsite)
;
IF KEYWORD_SET (prior) THEN BEGIN
FOR k=0,2 DO BEGIN
  openr,u,mdir0+'fpar'+string(k+1,fo='(I1)')+'_site.dat',/get_lun
  FOR s=0,nsite-1 DO BEGIN
    readf,u,lat,lon,nyr,year0
    FOR i=0,nyr-1 DO BEGIN
      readf,u,in & mfpar0[*,i,k,s]=in
    ENDFOR
  ENDFOR
  free_lun,u
ENDFOR
mfpar0=reform(mfpar0[*,yr0:yr1,*,*],(yr1-yr0+1)*365,maxpft,nsite)
END

IF KEYWORD_SET(ps) THEN BEGIN
  set_plot,'ps'
  device,file='fapar_sites.ps'
  IF ps EQ 2 THEN device, /color
ENDIF ELSE BEGIN
  ps = 0
ENDELSE

; TOTAL FAPAR COST FUNCTION
cftot = 0.
cftot0 = 0.

; go through subsites
FOR s = 0, nsite-1 DO BEGIN

; SATELLITE FAPAR
ncdf_read,d,file=dir+'fapar_redres_'+site[s]+'.nc',/all
x0=coord[0,s]-1 & x1=coord[1,s]-1
y0=coord[2,s]-1 & y1=coord[3,s]-1
fpar=total(total(d.fapar[x0:x1,y0:y1,*]/d.errfapar[x0:x1,y0:y1,*]^2,1),1)/$
  total(total(1./d.errfapar[x0:x1,y0:y1,*]^2,1),1)
errfpar=total(total(1/d.errfapar[x0:x1,y0:y1,*],1),1)/$
  total(total(1./d.errfapar[x0:x1,y0:y1,*]^2,1),1)
fpar[where(errfpar GT 1e3)]=!VALUES.F_NAN
; check length of FAPAR data set:
l=where(finite(fpar))
print, site[s]+': data from - to (DOY/YEAR) '+strcompress(fix(d.doy[l[0]]),/rem)+'/'+strcompress(fix(d.year[l[0]]),/rem)+'  '+strcompress(fix(d.doy[l[n_elements(l)-1]]),/rem)+'/'+strcompress(fix(d.year[l[n_elements(l)-1]]),/rem)
print, site[s]+', valid daily entries:',n_elements(l),'('+string(n_elements(l)*100./(l[n_elements(l)-1]-l[0]),fo='(F4.1,"%")')+')'
; COMPOSITE MODEL FAPAR AND COST FUNCTION
f=frac[*,s]
mfapar = mfpar[*,0,s]*f[0]+mfpar[*,1,s]*f[1]+mfpar[*,2,s]*f[2]
; FAPAR PART OF COST FUNCTION
l = where (finite(fpar))
cf = 0.5 * total ((fpar[l]-mfapar[l])^2 / errfpar[l]^2)
cftot = cftot + cf
; R.M.S. DEVIATION
rms = sqrt (total ((fpar[l]-mfapar[l])^2) / n_elements(l))

IF KEYWORD_SET (prior) THEN BEGIN
mfapar0 = mfpar0[*,0,s]*f[0]+mfpar0[*,1,s]*f[1]+mfpar0[*,2,s]*f[2]
cf0 = 0.5 * total ((fpar[l]-mfapar0[l])^2 / errfpar[l]^2)
cftot0 = cftot0 + cf0
rms0 = sqrt (total ((fpar[l]-mfapar0[l])^2) / n_elements(l))
print, 'cost function, r.m.s. error prior/post for ',site[s],': ', cf0, rms0, cf, rms
ENDIF ELSE $
print, 'cost function for ',site[s],': ', cf, rms

; PLOT COMPARISON
plot,fpar,xtitle='days since 1/1/2002',ytitle='FAPAR',title=site[s],yrange=[0,1],/nodata
IF ps EQ 2 THEN BEGIN
  errplot,fpar-errfpar,fpar+errfpar,thick=5.,color=180
ENDIF ELSE BEGIN
    errplot,fpar-errfpar,fpar+errfpar
ENDELSE
oplot,fpar,psym=7
oplot,mfapar
IF KEYWORD_SET (prior) THEN oplot,mfapar0,li=1

; PLOT INDIVIDUAL PFTs
l=where(f GT 0.)
IF KEYWORD_SET (pft_plot) AND l[0] NE -1 THEN BEGIN
FOR i = 0, n_elements(l)-1 DO BEGIN
  k = l[i]
  oplot,mfpar[*,k,s]*f[k],li=k+1
  plots,/normal,[0.2,0.3],[0.9-0.1*k,0.9-0.1*k],li=k+1
  xyouts,/normal,0.35,0.9-0.1*k,'PFT '+strcompress(pft[k,s],/rem)
ENDFOR
ENDIF

IF NOT KEYWORD_SET(ps) THEN WAIT, 5.0

ENDFOR

IF KEYWORD_SET(ps) THEN BEGIN
  device,/close
  set_plot,'x'
ENDIF

IF KEYWORD_SET (prior) THEN $
print, 'total cost function prior/post: ', cftot0, cftot $
ELSE $
print, 'total cost function : ', cftot

END

