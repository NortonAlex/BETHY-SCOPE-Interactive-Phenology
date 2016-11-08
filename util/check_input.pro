PRO check_input, name=name, lon = rlon, lat = rlat

path = '../input/climate/'

IF NOT KEYWORD_SET(name) THEN name = 'tmax'
IF NOT KEYWORD_SET(rlon) THEN rlon = 69.
IF NOT KEYWORD_SET(rlat) THEN rlat = 35.

ncdf_read,field,file=path+'forcing.hires.'+name+'.1979-2005.nc',variables=name,/no_struct

s = size (field) & nx = s[1] & ny = s[2] & nt = s[3]

lon = (findgen (nx) + 0.5) / nx * 360. - 180.
lat = (findgen (ny) + 0.5) / ny * 180. - 90.
t = (findgen(nt) + 0.5) / 365. + 1979

i = fix ((rlon + 180.) * nx / 360.)
j = fix ((rlat + 90.) * ny / 180.)

plot,t,field[i,j,*],/xstyle,title='lat: '+strcompress(rlat)+' lon: '+strcompress(rlon),ytitle=name

set_plot,'ps'
device,file=name+'_'+strcompress(rlat,/rem)+'_'+strcompress(rlon,/rem)+'.ps'
plot,t,field[i,j,*],/xstyle,title='lat: '+strcompress(rlat)+' lon: '+strcompress(rlon),ytitle=name
device,/close
set_plot,'x'

device,decompose=0
loadct,40
tvlct, 255B, 255B, 255B
miss = 1e20
mask = field[*,*,0] ne miss

STOP

END
