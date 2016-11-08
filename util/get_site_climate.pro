ncdf_read, file='input/climate/forcing.hires.tmax.1979-2005.nc',data,/all
tmax=data.tmax[*,*,7670:9130]
ncdf_read, file='input/climate/forcing.hires.tmin.1979-2005.nc',data,/all
tmin=data.tmin[*,*,7670:9130]
ncdf_read, file='input/climate/forcing.hires.precip.1979-2005.nc',data,/all
precip=data.precip[*,*,7670:9130]
ncdf_read, file='input/climate/forcing.hires.swdown.1979-2005.nc',data,/all
swdown=data.swdown[*,*,7670:9130]

openw,1,'forcing.sites.2000-2003.txt'

;loobos
x=93
y=71
printf,1,'loobos'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

;sodankyla
x=103
y=78
printf,1,'sodankyla'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

;zotino
x=134
y=75
printf,1,'zotino'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

;hainich
x=95
y=70
printf,1,'hainich'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

;Manaus
x=59
y=43
printf,1,'manaus'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

;Tapajos
x=62
y=43
printf,1,'tapajos'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

;Maun
x=101
y=35
printf,1,'maun'
printf,	1,tmax[x,y,*],form='(1461f8.2)'
printf,	1,tmin[x,y,*],form='(1461f8.2)'
printf,	1,precip[x,y,*],form='(1461f8.2)'
printf,	1,swdown[x,y,*],form='(1461f8.2)'

close,1

.run
n=0
for i=1979,2006 do begin
print,i,n
if ((i mod 4) eq 0) then begin
n=n+366
endif else begin
n=n+365
endelse
endfor
end


