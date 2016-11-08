prior=fltarr(1464)
filen='cost_fapar_prior.dat'
openr, 1, filen
readf, 1, prior, format='(366f6.3)'
close, 1

post=fltarr(1464)
filen='cost_fapar_post.dat'
openr, 1, filen
readf, 1, post, format='(366f6.3)'
close, 1

obs=fltarr(1464)
filen='cost_fapar_obs.dat'
openr, 1, filen
readf, 1, obs, format='(366f6.3)'
close, 1

nd = 4 * 366
yr0 = 2000
x = dindgen(nd) + 1.
x = yr0 + (x - 0.5)/366.

c=getcolor(/load)

obsv=obs[where(obs gt 0.)]
xv=x[where(obs gt 0.)]

plot,xv,obsv,xrange=[2000,2004],yrange=[0,1],psym=1
errplot,xv,obsv-.1,obsv+.1
oplot,x,post,color=c.blue
oplot,x,prior,color=c.red


set_plot,'ps'
device,file='loobos_fapar_optimisation_2pft.ps',/color
plot,xv,obsv,xrange=[2000,2004],yrange=[0,1],psym=1,xtitle='Year',ytitle='FAPAR'
errplot,xv,obsv-.1,obsv+.1
oplot,x,post,color=c.blue
oplot,x,prior,color=c.red
device,/close
set_plot,'x'

