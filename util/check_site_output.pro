openr, fu, 'fpar_site.dat', /get_lun
readf, fu, lat, lon, nyr, yr0
nd = 12*nyr
fpar = fltarr(nd)
readf, fu, fpar              
close,fu

x = findgen(nd) + 1.
x = yr0 + (x - 0.5)/12.

plot, x, fpar, xtitle = 'Year', ytitle = 'Fpar'

set_plot,'ps'
device,file='hainich_fpar.ps'
plot,x, fpar, xtitle = 'Year', ytitle = 'Fpar'
device,/close
set_plot,'x'


openr, fu, 'lai_site.dat', /get_lun
readf, fu, lat, lon, nyr, yr0
lai = fltarr(nd)
readf, fu, lai
close,fu

plot, x, lai, xtitle = 'Year', ytitle = 'LAI [m2/m2]'

set_plot,'ps'
device,file='hainich_lai.ps'
plot,x, lai, xtitle = 'Year', ytitle = 'LAI [m2/m2]'
device,/close
set_plot,'x'


openr, fu, 'npp_site.dat', /get_lun
readf, fu, lat, lon, nyr, yr0
npp = fltarr(nd)
readf, fu, npp
close,fu

plot, x, npp, xtitle = 'Year', ytitle = 'NPP [gC/m2]'

set_plot,'ps'
device,file='hainich_npp.ps'
plot,x, npp, xtitle = 'Year', ytitle = 'NPP [gC/m2]'
device,/close
set_plot,'x'
