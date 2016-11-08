PRO bethy_lores_grid

; Version 2, 26 June 2007
; correction on 27 July 2007
; Version 3, 23 October 2003: changed albedo for areas identified as 'ice'; now on 'ppm' sandbox
; Version 3.1, 26 October 2007: correct points with with total PFT fraction equal to 0


; IDL programme to generate input file for CCDAS in TM2 resolution (="lo res")
; it reads other information from a earlier CCDAS input file using only
; IMBETHY and generates a new file adding various new field
; identical to bethy_grid.pro, but with additional conversion of bucket
; size to "lores" (also, 'mtype' does not exist for the lores input file)

fullBETHYfile = '/Users/wolfgang/Models/fullBETHY/output/bethy_VICA/vegtypes_VIC_2deg.dat'
fullBETHYformat = '(f6.1,f7.1,5i3,3f7.4,3f8.2,3f6.3,3f6.2,f6.2,f8.2,f6.3)'
IMBETHY_infile = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/adbethy_loresgrid3veg_79-00.nc'
newCCDAS_outfile = '/Users/wolfgang/Models/CCDAS/ppm/control_bethy/bethy_loresgrid3veg_v2.nc'
surface_file = '/Users/wolfgang/Models/fullBETHY/input/surface6_VIC_2deg.dat'
surface_format = '(f6.1, f7.1, i5, f7.1, i4, i3, i3, i3, i3, f6.3, f7.1)'
FAOlayers_file = '/Users/wolfgang/Models/fullBETHY/input/avgprofiles.large'
missing = -9999.

; get 'wmax' (maximum plant-available soil moisture)
loadm, fullBETHYfile, v, format=fullBETHYformat, gridwidth=2.0, /center, co=20, oc=-1e6
wmaxmin = 10.
wmax = v[*,*,18]
land = wmax NE -1e6
wmax = (wmax > wmaxmin)*land - 1e6 * (land EQ 0)
ncdf_read, file=IMBETHY_infile, /all, in

; get FAO soil type, surface reflectance when dry or wet after Wilson & Henderson-Sellers
loadm, surface_file, s, format = surface_format, gridwidth=2.0, /center, co=9, ocean=-1e6
land = s[*,*,0] ne -1e6
FAO = s[*,*,2] > 0
st = s[*,*,3] > 0
ice = st EQ 34
water = st EQ 0
; adopted from BETHY source code 'soil.f', subroutine 'getsoildata' and 'bethy.options'
rsoil = [0.070,  0.150,  0.100,  0.200,  0.180,  0.350]
; default is light soilisb = st*0 + 5
isb (where ((st GE 17 AND st LE 22) OR st EQ 30)) = 3
isb (where ((st GE 23 AND st LE 28) OR st EQ 31)) = 1
rhosw = rsoil(isb-1)
rhosd = rsoil(isb)
;rhosw(where(ice)) = 0.70
;rhosd(where(ice)) = 0.70
rhosw(where(ice)) = rsoil[2] ; ice areas set to medium bright soil
rhosd(where(ice)) = rsoil[3]
rhosw(where(water)) = 0.02
rhosd(where(water)) = 0.02
rhosw = rhosw * land - 1e6 * (land eq 0)
rhosd = rhosd * land - 1e6 * (land eq 0)
; determine parameters for soil evaporation model, adopted from 'exsoildata.f'
nprof = 134
maxhor = 4
nhor=intarr(nprof)
scode = strarr(nprof)
sand=fltarr(maxhor,nprof) & silt=fltarr(maxhor,nprof) & clay=fltarr(maxhor,nprof)
thick=fltarr(maxhor,nprof)
line = ' '
openr, 1, FAOlayers_file
for ip = 0, nprof-1 do begin
  readf, 1, format = '(A2,I2)', line, nhor0
  nhor[ip] = nhor0
  scode[ip] = line
  for ih = 0, maxhor-1 do begin
    readf, 1, format = '(9X,4F8.4)', thick0, sand0, silt0, clay0
    thick[ih,ip]=thick0
    sand[ih,ip]=sand0 & silt[ih,ip]=silt0 & clay[ih,ip]=clay0
  endfor
  readf, 1, line
  print, line
endfor
close, 1
sand1 = reform (sand[0,FAO], n_elements(FAO[*,0]), n_elements(FAO[0,*])) * land
sand1(where(sand1 eq -1.)) = 0.8 ; generic texture, see 'exsoildata.f'
des = (5.62 - 2.56 * sand1) * land - 1e6 * (land eq 0)
evap1 = (14.29 - 9.23 * sand1) * land - 1e6 * (land eq 0)

; correct fractional cover where the sum is 0
l = where(total(in.frac,3) eq 0)
if l[0] ne -1 then in.frac[l] = 1. ; set first entry to 100% of grid cell

; transform fields to lores
missing = -9999.
lsm_lores = in.gridnum ne missing
wmax_lores = bethy2tm (wmax, miss=-1e6) * lsm_lores + missing * (lsm_lores eq 0)
rhosw_lores = bethy2tm (rhosw, miss=-1e6) * lsm_lores + missing * (lsm_lores eq 0)
rhosd_lores = bethy2tm (rhosd, miss=-1e6) * lsm_lores + missing * (lsm_lores eq 0)
des_lores = bethy2tm (des, miss=-1e6) * lsm_lores + missing * (lsm_lores eq 0)
evap1_lores = bethy2tm (evap1, miss=-1e6) * lsm_lores + missing * (lsm_lores eq 0)

; check fields visually
device,decompose=0
loadct,40
tvlct, 255B, 255B, 255B
;tvf,FAO,/bar,/map,title='FAO soil type'
;tvf,sand1*100,/bar,/map,title='% sand at surface'
window,0 & tvf,wmax_lores,/bar,/map,missing=missing,title='max plant avail soil water [mm]'
window,1 & tvf,rhosw_lores,/bar,missing=missing,/map,title='wet soil albedo'
window,2 & tvf,rhosd_lores,/bar,missing=missing,/map,title='dry soil albedo'
window,3 & tvf,des_lores,/bar,missing=missing,/map,title='desorptivity [mm/sqrt(day)]'
window,4 & tvf,evap1_lores,/bar,missing=missing,/map,title='phase 1 evaporation [mm]'

; create NetCDF file

id = NCDF_CREATE(newCCDAS_outfile, /CLOBBER)

  ; define dimensions
dlon   = NCDF_DIMDEF(id, 'lon', n_elements(in.lon))
dlat   = NCDF_DIMDEF(id, 'lat', n_elements(in.lat))
dvnum  = NCDF_DIMDEF(id, 'vnum', 3)
dtime  = NCDF_DIMDEF(id, 'time', /UNLIMITED)


vlon = NCDF_VARDEF(id, 'lon', [dlon], /FLOAT)
NCDF_ATTPUT, id, vlon, 'long_name', 'longitude', /CHAR
NCDF_ATTPUT, id, vlon, 'units', 'degrees_east', /CHAR

vlat = NCDF_VARDEF(id, 'lat', [dlat], /FLOAT)
NCDF_ATTPUT, id, vlat, 'long_name', 'latitude', /CHAR
NCDF_ATTPUT, id, vlat, 'units', 'degrees_north', /CHAR

vvnum =  NCDF_VARDEF(id, 'vnum', [dvnum], /LONG)
NCDF_ATTPUT, id, vvnum, 'long_name', 'vegetation number', /CHAR
NCDF_ATTPUT, id, vvnum, 'units', '', /CHAR

vtime = NCDF_VARDEF(id, 'time', [dtime], /FLOAT)
NCDF_ATTPUT, id, vtime, 'long_name', 'time', /CHAR
NCDF_ATTPUT, id, vtime, 'units', 'month', /CHAR

   ; define variables
  
vgridnum = NCDF_VARDEF(id, 'gridnum', [dlon, dlat], /LONG)
NCDF_ATTPUT, id, vgridnum, 'long_name', 'land grid cell number starting at 180E/90N', /CHAR
NCDF_ATTPUT, id, vgridnum, 'units', '', /CHAR
NCDF_ATTPUT, id, vgridnum, 'missing_value', '-9999', /LONG

velev = NCDF_VARDEF(id, 'elev', [dlon, dlat], /LONG)
NCDF_ATTPUT, id, velev, 'long_name', 'Mean grid cell elevation', /CHAR
NCDF_ATTPUT, id, velev, 'units', 'm', /CHAR
NCDF_ATTPUT, id, velev, 'missing_value', '-9999', /LONG

vtype = NCDF_VARDEF(id, 'type', [dlon, dlat, dvnum], /LONG)
NCDF_ATTPUT, id, vtype, 'long_name', 'Vegetation type (0-13)', /CHAR
NCDF_ATTPUT, id, vtype, 'units', '', /CHAR
NCDF_ATTPUT, id, vtype, 'missing_value', '-9999', /LONG

vfrac = NCDF_VARDEF(id, 'frac', [dlon, dlat, dvnum], /FLOAT)
NCDF_ATTPUT, id, vfrac, 'long_name', 'Vegetation type fraction', /CHAR
NCDF_ATTPUT, id, vfrac, 'units', '', /CHAR
NCDF_ATTPUT, id, vfrac, 'missing_value', '-9999', /FLOAT

vwmax = NCDF_VARDEF(id, 'wmax',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vwmax, 'long_name', 'Maximum plant-available soil moisture', /CHAR
NCDF_ATTPUT, id, vwmax, 'units', 'mm', /CHAR
NCDF_ATTPUT, id, vwmax, 'missing_value', '-9999', /FLOAT

vrhosw = NCDF_VARDEF(id, 'rhosw',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vrhosw, 'long_name', 'Wet-soil albedo', /CHAR
NCDF_ATTPUT, id, vrhosw, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vrhosw, 'missing_value', '-9999', /FLOAT

vrhosd = NCDF_VARDEF(id, 'rhosd',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vrhosd, 'long_name', 'Dry-soil albedo', /CHAR
NCDF_ATTPUT, id, vrhosd, 'units', 'fraction', /CHAR
NCDF_ATTPUT, id, vrhosd, 'missing_value', '-9999', /FLOAT

vdes = NCDF_VARDEF(id, 'des',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vdes, 'long_name', 'Soil water desortptivity', /CHAR
NCDF_ATTPUT, id, vdes, 'units', 'mm/sqrt(day)', /CHAR
NCDF_ATTPUT, id, vdes, 'missing_value', '-9999', /FLOAT

vevap1 = NCDF_VARDEF(id, 'evap1',[dlon,dlat], /FLOAT)
NCDF_ATTPUT, id, vevap1, 'long_name', 'Phase 1 soil evaporation', /CHAR
NCDF_ATTPUT, id, vevap1, 'units', 'mm', /CHAR
NCDF_ATTPUT, id, vevap1, 'missing_value', '-9999', /FLOAT

vmpot = NCDF_VARDEF(id, 'mpot', [dlon, dlat, dtime], /FLOAT)
NCDF_ATTPUT, id, vmpot, 'long_name', 'Mean potential shortwave radiation', /CHAR
NCDF_ATTPUT, id, vmpot, 'units', 'W/m2', /CHAR
NCDF_ATTPUT, id, vmpot, 'missing_value', '-9999', /FLOAT

   ; define global attributes
NCDF_ATTPUT, id, /GLOBAL, 'title', 'Boundary conditions for CCDAS-BETHY', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'author', 'Wolfgang Knorr & Marko Scholze', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'affiliation', 'QUEST, Uni Bristol', /CHAR
NCDF_ATTPUT, id, /GLOBAL, 'date', 'Version 2, creation date 26 June 2007', /CHAR

   ; change to DATA mode and put data into file
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, vlon, in.lon
NCDF_VARPUT, id, vlat, in.lat
NCDF_VARPUT, id, vvnum, indgen(3)+1
NCDF_VARPUT, id, vtime, indgen(n_elements(in.mpot[0, 0, *]))+1
NCDF_VARPUT, id, vgridnum, fix(in.gridnum)
NCDF_VARPUT, id, velev, in.elev
NCDF_VARPUT, id, vtype, in.type
NCDF_VARPUT, id, vfrac, in.frac
NCDF_VARPUT, id, vwmax, wmax_lores
NCDF_VARPUT, id, vrhosw, rhosw_lores
NCDF_VARPUT, id, vrhosd, rhosd_lores
NCDF_VARPUT, id, vdes, des_lores
NCDF_VARPUT, id, vevap1, evap1_lores
NCDF_VARPUT, id, vmpot, in.mpot

   ; close file
NCDF_CLOSE, id

END


