PRO add_site_data, infile

IF NOT KEYWORD_SET (infile) THEN infile = 'site_specs_multiple'

path = '../control_bethy/'
line = ' '
name = ' '
vtype = intarr(3)
vfrac = fltarr(3)
coord = intarr(4)

ncdf_read,file=path+'bethy_grid3veg_v2.nc',/all,grid
nx = n_elements(grid.lon)
ny = n_elements(grid.lat)

openr,u,path+infile+'.dat',/get_lun
openw,v,path+infile+'_new.dat',/get_lun

readf,u,line
printf,v,'lat lon elev pft frac scene_x_from scene_x_to scene_y_from scene_y_to wmax rhosw rhosd des evap1 site_name'

WHILE NOT EOF(u) DO BEGIN

readf,u,format='(2f7.2,i6,3i3,3f5.2,4I4,1X,A)', lat, lon, elev, vtype, vfrac, coord, name

lons = grid.lon#(fltarr(ny)+1)*(grid.gridnum ne -9999)-9999*(grid.gridnum eq -9999)
lats = (fltarr(nx)+1)#grid.lat*(grid.gridnum ne -9999)-9999*(grid.gridnum eq -9999)
l = sort((lon-lons)^2+(lat-lats)^2) & l=l[0] ; find nearest land cell
x = l mod nx & y = l / nx
print,'lon: ',lon,lons[l],grid.lon[x],'      lat: ',lat,lats[l],grid.lat[y]
;x = sort(abs(lon-grid.lon)) & x = x[0]
;y = sort(abs(lat-grid.lat)) & y = y[0]
;print,'lon: ',lon,grid.lon[x],'      lat: ',lat,grid.lat[y]

printf,v,format='(2f7.2,i6,3i3,3f5.2,4I4,F7.1,2F5.2,F6.2,F7.1,1X,A)', lat, lon, elev, vtype, vfrac, coord, grid.wmax[x,y], grid.rhosw[x,y], grid.rhosd[x,y], grid.des[x,y], grid.evap1[x,y], name

ENDWHILE

free_lun, u

END

