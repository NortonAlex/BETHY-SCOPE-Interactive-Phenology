PRO check_input_revNov2008, name=name, smin = smin

path = '/Users/wolfgang/Models/VIC_drivers/revision_nov_2008/'
missing = 1e20
year0 = 1979
year1 = 2008
device,decompose=0
loadct,40
tvlct, 255B, 255B, 255B

IF NOT KEYWORD_SET(name) THEN name = 'tmax'
IF NOT KEYWORD_SET(smin) THEN smin = 20 ; standard deviation threshold

ncdf_read,field,file=path+'forcing.2x2.'+name+'.'+strcompress(year0,/rem)+'-'+strcompress(year1,/rem)+'.nc',variables=name,/no_struct
s = size (field) & nx = s[1] & ny = s[2] & nt = s[3]

lon = (findgen (nx) + 0.5) / nx * 360. - 180.
lat = (findgen (ny) + 0.5) / ny * 180. - 90.
t = (findgen(nt) + 0.5) / 365. + 1979
land = field[*,*,0] ne missing

; compute monthly climatology of mean and standard deviation
dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
nyr = year1 - year0 + 1
yr = year0
mth = 0
dom = 0 ; day-of-month
mfield = fltarr (nx, ny, 12)
pmth = intarr (nt) ; pointer to month
pyr = intarr (nt) ; pointer to year
for day = 0, nt - 1 do begin
  dom = dom + 1
  dpmth = dpm[mth]
  if mth eq 1 and yr mod 4 eq 0 then dpmth = 29 ; leap year
  mfield[*,*,mth] = mfield[*,*,mth] + field[*,*,day] / dpmth / nyr
  pmth[day] = mth & pyr[day] = yr-year0
  ;print, 'day yr mth: ', day, yr, mth
  if dom GE dpmth then begin
    dom = 0
    mth = mth + 1
  endif
  if mth EQ 12 then begin
    mth = 0
    yr = yr + 1
  endif
endfor
for mth = 0, 11 do mfield[*,*,mth] = mfield[*,*,mth] * land + missing * (land eq 0)

; subtract climatology from field (seasonal adjustment)
for day = 0, nt-1 do field[*,*,day] = field[*,*,day] - mfield[*,*,pmth[day]]

; compute 3-D field of standard deviation of seasonally adjusted field
sfield = fltarr (nx, ny)
nfield = intarr (nx, ny)
for x = 0, nx-1 do begin
  for y = 0, ny-1 do begin
    if field[x,y,0] ne missing then begin
      sfield[x,y] = stddev(field[x,y,*])
      nfield[x,y] = total (abs(field[x,y,*]/sfield[x,y]) gt smin)
    endif
  endfor
  print, x
endfor
nfield = nfield + missing * (land eq 0)

window

tvf,sfield+missing*(land eq 0),miss=missing,/bar,title=name+': std. dev. of daily value after monthly anomaly correction',tif=name+'_astd.tif'

tvf,nfield,miss=missing,/bar,title=name+': number of days with abs. anomaly > '+strcompress(smin,/rem)+' std. dev. of anomaly time series',tif=name+'_'+strcompress(smin,/rem)+'s.tif'
; number of > smin sigma deviating points per month and year

ndev = intarr (12, nyr)
for day = 0, nt-1 do ndev[pmth[day], pyr[day]] = ndev[pmth[day], pyr[day]] + total (abs(field[*,*,day]/(sfield[*,*]+(sfield[*,*] eq 0))) gt smin)

set_plot,'ps'
device,file=name+'_'+strcompress(smin,/rem)+'s.ps'
plot,findgen(12*nyr)/12+year0,reform(ndev,12*nyr),/xstyle,title='cases with '+strcompress(name,/rem)+' > '+strcompress(smin,/rem)+' stddev from monthly climatology (per month)'
device,/close
set_plot,'x'

l = where(nfield*land gt 0)
x = l mod nx
y = fix (l / nx)

IF l[0] NE -1 THEN BEGIN

set_plot,'ps'
device,file=name+'_points_'+strcompress(smin,/rem)+'s.ps'
FOR i = 0, n_elements(l)-1 DO $
  plot,field[x[i],y[i],*]+mfield[x[i],y[i],pmth],title=name+' daily value: '+strcompress(lat[y[i]],/rem)+' lat, '+strcompress(lon[x[i]],/rem)+' lon'
device,/close
set_plot,'x'

ENDIF

;STOP

END
