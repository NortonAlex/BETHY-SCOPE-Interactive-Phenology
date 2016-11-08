openr,1,'output_prior/nwr_modconc.dat'
read_graph,1,h,x,yp
close,1

openr,1,'output/nwr_modconc.dat'
read_graph,1,h,x,y
close,1


print,h
stat,y[0,*]
stat,y[1,*]
b = 12
filter_data,x=reform(x[b:*]), y=reform( y[0,b:*]), cutoff1=80., cutoff2=650.,$
			 npoly=3, nharm=4,  interval=7., tzero=1979.,$
			 /even, climseas= obs

filter_data,x=reform(x[b:*]), y=reform( yp[1,b:*]), cutoff1=80., cutoff2=650.,$
			 npoly=3, nharm=4,  interval=7., tzero=1979.,$
			 /even, climseas= prior

filter_data,x=reform(x[b:*]), y=reform( y[1,b:*]), cutoff1=80., cutoff2=650.,$
			 npoly=3, nharm=4,  interval=7., tzero=1979.,$
			 /even, climseas= model


stat,obs[1,*]
stat,model[1,*]
stat,model[1,*]-mean(y[2,*] < 100)
stat,model[1,*]+mean(y[2,*] < 100)
meanunc = mean(y[2, *] * (y[2,*] lt 100))
u = ['Jan', 'Feb', 'Mar', 'Apr', 'May','Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
c = getcolor(/load)

!p.font=0
set_plot, 'ps'
device, filename='mlo_seas.eps', encapsulated=0,xoffset=2.,yoffset=2.0,xsize=8.5,ysize=6.1, SET_FONT = 'Helvetica',font_size=8, /color

plot,obs[0,*],obs[1,*],xrange=[1990.,1991.],yrange=[-10.,10.],thick=2, xticks = 11, xstyle = 1, xtickn = u, xtitle = 'Month', ytitle = '[ppm]'
oplot,obs[0,*],obs[1,*]-meanunc,lines=1
oplot,obs[0,*],obs[1,*]+meanunc,lines=1
oplot,obs[0,*],model[1,*],th = 2, color = c.red
oplot,obs[0,*],prior[1,*],th = 2, color = c.blue

device, /close
set_plot, 'x'
