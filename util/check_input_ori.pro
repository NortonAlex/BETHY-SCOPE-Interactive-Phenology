PRO check_input_ori, name=name, lon = rlon, lat = rlat

path = '~/Data/VIC/'
missing = 1e20
year0 = 1979
device,decompose=0
loadct,40
tvlct, 255B, 255B, 255B

IF NOT KEYWORD_SET(name) THEN name = 'tmax'
IF NOT KEYWORD_SET(rlon) THEN rlon = 69.
IF NOT KEYWORD_SET(rlat) THEN rlat = 35.

ncdf_read,field,file=path+'forcing.2x2.'+name+'.1979-2005.nc',variables=name,/no_struct

s = size (field) & nx = s[1] & ny = s[2] & nt = s[3]

lon = (findgen (nx) + 0.5) / nx * 360. - 180.
lat = (findgen (ny) + 0.5) / ny * 180. - 90.
t = (findgen(nt) + 0.5) / 365. + 1979
land = field[*,*,0] ne missing

i = fix ((rlon + 180.) * nx / 360.)
j = fix ((rlat + 90.) * ny / 180.)

; compute monthly climatology of mean and standard deviation
dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
nyr = nt / 365
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
  ;print, day, yr, mth
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

smin = 6 ; standard deviation threshold
smin = 10 ; standard deviation threshold

; compute 3-D field of standard deviation of seasonally adjusted field
sfield = fltarr (nx, ny)
nfield = intarr (nx, ny)
for x = 0, nx-1 do begin
  for y = 0, ny-1 do begin
    if field[x,y,0] ne missing then begin
      sfield[x,y] = stddev(field[x,y,*])
      nfield[x,y] = total (field[x,y,*]/sfield[x,y] gt smin)
    endif
  endfor
  print, x
endfor

; number of > smin sigma deviating points per month and year
ndev = intarr (12, nyr)
for day = 0, nt-1 do ndev[pmth[day], pyr[day]] = ndev[pmth[day], pyr[day]] + total (abs(field[*,*,day]/(sfield[*,*]+(sfield[*,*] eq 0))) gt smin)

plot,findgen(12*nyr)/12+year0,reform(ndev,12*nyr),/xstyle,title='cases with '+strcompress(name,/rem)+' > '+strcompress(smin,/rem)+' stddev from monthly climatology (per month)'

;tvf,field[*,*,day]/(sfield[*,*]+(sfield[*,*] eq 0)),/bar

; deviation from climatology
plot,t,field[i,j,*],/xstyle,title='lat: '+strcompress(rlat)+' lon: '+strcompress(rlon),ytitle=name

set_plot,'ps'
device,file=name+'_dev_'+strcompress(rlat,/rem)+'_'+strcompress(rlon,/rem)+'.ps'
plot,t,field[i,j,*],/xstyle,title='lat: '+strcompress(rlat)+' lon: '+strcompress(rlon)+' dev. from monthly clim.',ytitle=name
device,/close
set_plot,'x'

set_plot,'ps'
device,file=name+'_count_'+strcompress(smin,/rem)+'sigma.ps'
plot,findgen(12*nyr)/12+year0,reform(ndev,12*nyr),/xstyle,title='cases with '+strcompress(name,/rem)+' > '+strcompress(smin,/rem)+' stddev from monthly climatology (per month)'
device,/close
set_plot,'x'


STOP

END
