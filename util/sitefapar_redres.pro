PRO sitefapar_redres

; routine to write site FAPAR data to NetCDF
; to retrieve the data, use:
; ncdf_read,file='fapar_redres_[site name].nc',data,/all
; help,/str,data
; site_name=string(data.name)                     
; help,site_name
; (using Martin Schultz's ncdf_read.pro)

site_list = ['Botswana', 'Manaus', 'Sodankyla', 'Loobos', 'Tapajos', 'Zotino', 'Hainich']
output_site_list = ['BW-Ma1', 'manaus', 'FI-Sod', 'NL-Loo', 'tapajos', 'RU-Zot', 'hainich']

file='Daily.dat'
indir = '/Users/wolfgang/Documents/Projects/RS-CCDAS/data/low resolution site scenes/'
outfile = 'fapar_redres'
outdir = '/Users/wolfgang/Models/CCDAS/input/eddy_sites/'
;
nsite = n_elements(site_list)
scale_fapar =0.003937007859349251
offset_fapar =-0.003937007859349251
fapar_missing = 0.0
errfapar_standard = 0.1
errfapar_missing = 1e9
;
startyear = 2002
endyear = 2003
daysperyear = 365
nyear = endyear-startyear+1
ntime = daysperyear*nyear
nlat = 15
nlon = 15
fapar = fltarr (nlon,nlat,ntime,nsite) + fapar_missing
errfapar = fltarr (nlon,nlat,ntime,nsite) + errfapar_missing
;
FOR site = 0, nsite-1 DO BEGIN
  namefile=indir+'ReducedResolution/'+site_list[site]+'/'+file
  data = LoadProfileData_RedRes(namefile)
  FOR i = 0, n_elements(data)-1 DO BEGIN
    time = data[i].startday + (data[i].year-startyear)*daysperyear
    fapar(*,*,time,site) = offset_fapar + scale_fapar * data[i].fapar
;    errfapar(*,*,time,site) = errfapar_standard * (data[i].flag EQ 101B) $
;                            + errfapar_missing * (data[i].flag NE 101B)
    errfapar(*,*,time,site) = errfapar_standard * $
                              (data[i].fapar GT 0B AND (data[i].flag EQ 101B OR data[i].flag EQ 102B)) $
                            + errfapar_missing * $
                              (data[i].fapar EQ 0B OR (data[i].flag NE 101B AND data[i].flag NE 102B))
  ENDFOR
  print, 'Site ',site_list[site],' processed.'
ENDFOR
doy = intarr(ntime)
year = intarr(ntime)
for i=0,nyear-1 do doy[i*daysperyear:(i+1)*daysperyear-1]=indgen(daysperyear)+1
for i=0,nyear-1 do year[i*daysperyear:(i+1)*daysperyear-1]=startyear+i
;
loadct,41
device,decompose=0
window,0,xs=500,ysize=400
window,1,xs=500,ysize=400
site = nsite-1
for i=0,ntime-1 do begin
  wset, 0
  tvf,fapar[*,*,i,site],/rev,/bar,title='FAPAR   day:'+strcompress(i)+'   SITE:'+strcompress(site_list[site]),min=0,max=0.8
  wset, 1
  tvf,errfapar[*,*,i,site],/rev,/bar,title='FAPAR   day:'+strcompress(i)+'   SITE:'+strcompress(site_list[site]),min=0,max=0.8
  wait, 0.2
end
; FILTER OUTLIERS USING RUNNING MEAN
; cannot use 'smooth', as unable to deal with non-valid points
; instead, compute average with uncertainties for 7 entries around current
; current value, then subtract average from current value and divide by
; combined uncertainty. If this values is below -3,
; treat as cloud contaminated outlier and remove.
f=fapar/errfapar^2 ; weighted fapar
e=1./errfapar ; weighted error
w=1./errfapar^2    ; the weights
;af = fltarr(nlon,nlat,ntime,nsite)
;ef = fltarr(nlon,nlat,ntime,nsite)
df = fltarr(nlon,nlat,ntime,nsite) ; array to store the normalized deviation
nav = 7 ; averaging interval (must be odd integral)
threshold = -2 ; threshold negative outlier in units of uncertainty
FOR site = 0, nsite - 1 DO BEGIN
  FOR i = 0, nlon - 1 DO BEGIN
    FOR j = 0, nlat - 1 DO BEGIN
      FOR t = (nav-1)/2, ntime-1-(nav-1)/2 DO BEGIN
        t0 = t - (nav-1)/2
        t1 = t + (nav-1)/2
  ;af[i,j,t,site] = total (f[i,j,t0:t1,site]) / total (w[i,j,t0:t1,site])
  ;ef[i,j,t,site] = total (e[i,j,t0:t1,site]) / total (w[i,j,t0:t1,site])
  ;df[i,j,t,site] = (fapar[i,j,t,site] - af[i,j,t,site]) / ef[i,j,t,site]
        df[i,j,t,site] = (fapar[i,j,t,site] - $
                    total (f[i,j,t0:t1,site]) / total (w[i,j,t0:t1,site])) $
                 / (total (e[i,j,t0:t1,site]) / total (w[i,j,t0:t1,site]))
      ENDFOR
    ENDFOR
  ENDFOR
ENDFOR
; mask out these points:
mask = (df LE threshold AND errfapar LE 0.1)
errfapar = errfapar * (mask EQ 0) + errfapar_missing * mask
print,'Data filter applied. Total points removed: ',fix(total(mask))
;
; create NetCDF files
FOR site = 0, nsite-1 DO BEGIN
;
filename=outdir+outfile+'_'+output_site_list[site]+'.nc'
id = NCDF_CREATE(filename, /CLOBBER)
; define dimensions
dlon  = NCDF_DIMDEF(id,'lon',nlon)
dlat  = NCDF_DIMDEF(id,'lat',nlat)
dtime = NCDF_DIMDEF(id,'time',ntime)
; define dimensions
vdoy = NCDF_VARDEF(id, 'DOY', [dtime], /FLOAT)
NCDF_ATTPUT, id, vdoy, 'long_name', 'day of year', /CHAR
NCDF_ATTPUT, id, vdoy, 'units', 'days', /CHAR
vyear = NCDF_VARDEF(id, 'year', [dtime], /FLOAT)
NCDF_ATTPUT, id, vyear, 'long_name', 'calendar year', /CHAR
NCDF_ATTPUT, id, vyear, 'units', 'years', /CHAR
; define variables
vfapar = NCDF_VARDEF(id, 'fapar', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vfapar, 'long_name', 'fraction of absorbed photosynhetically active radation (FAPAR)', /CHAR
NCDF_ATTPUT, id, vfapar, 'units', '', /CHAR
verrfapar = NCDF_VARDEF(id, 'errfapar', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, verrfapar, 'long_name', 'error in FAPAR', /CHAR
NCDF_ATTPUT, id, verrfapar, 'units', '', /CHAR
; define global attributes
NCDF_ATTPUT, id, /GLOBAL, 'title', 'Site Data (reduced resolution) for RS-CCDAS', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'site', output_site_list[site], /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'date', 'Version 2.0, creation date 17 Jan 2008', /CHAR
; change to DATA mode and put data into file
NCDF_CONTROL, id, /ENDEF
NCDF_VARPUT, id, vdoy, doy
NCDF_VARPUT, id, vyear, year
NCDF_VARPUT, id, vfapar, fapar[*,*,*,site]
NCDF_VARPUT, id, verrfapar, errfapar[*,*,*,site]
; close file
NCDF_CLOSE, id
print,'FAPAR file generated: ',filename
;
ENDFOR

END
