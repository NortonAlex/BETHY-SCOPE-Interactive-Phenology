openr,1,'h.dat'

opt = fltarr(3, 32)
readf, 1, opt
close, 1

c=getcolor(/load)

plot,opt[0,*],opt[1,*],yrange=[2200,2800],ystyle=4, xtitle='Iteration'
plot,opt[0,*],opt[2,*],yrange=[-1,200],/noerase,ystyle=9,ytitle='Gradient'

plot,opt[0,*],opt[1,*],yrange=[2200,2800], ytitle='Costfunction'


set_plot,'ps'
device,file='costf_grad.ps',/color
plot,opt[0,*],opt[1,*],yrange=[2200,2800],ystyle=4, xtitle='Iteration',/nodata
oplot,opt[0,*],opt[1,*],color=c.red
plot,opt[0,*],opt[2,*],yrange=[-1,200],/noerase,ystyle=9,ytitle='Gradient'
device,/close
set_plot,'x'


set_plot,'ps'
device,file='costf_grad_axis.ps',/color
plot,opt[0,*],opt[1,*],yrange=[2200,2800], ytitle='Costfunction', xtitle='Iteration'
device,/close
set_plot,'x'
