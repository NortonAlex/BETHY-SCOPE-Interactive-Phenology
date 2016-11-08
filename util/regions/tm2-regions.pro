ncdf_read, grid, filename = '/net/eiche/scratch/local2/CCDAS/ADBETHY/input/adbethy_loresgrid3veg_79-00.nc', /all
makearea, a, lon = 10., lat = 39.13
ar = fltarr(36, 24)
ar[*, 0] = a[*, 0]
ar[*, 23] = a[*, 45]
ar[*, 1:22] = rebin3(a[*, 1:44], 22, 2, /keepsum)
land = grid.elev ge 0.
ar = ar * land 

lon = (findgen(36)+ 0.5) / 36 * 360. - 180.
lat = (findgen(24)+ 0.5) / 24 * 180. - 90.

nmask=6
mask=fltarr(36,24,nmask)
; N-America
mask(*,*,0) = (lon GE -180 AND lon LE -20) # (lat GE 15 AND lat LE 90) *ar
; S-America
mask(*,*,1) = (lon GE -180 AND lon LE -20) # (lat GE -90 AND lat LE 15) * ar
; Europe
mask(*,*,2) = (lon GE -10 AND lon LE 45) # (lat GE 40 AND lat LE 90) * ar
; Asia
mask(*,*,3) = (lon GE 45 AND lon LE 180) # (lat GE 15 AND lat LE 90) * ar
; Africa
mask(*,*,4) = (lon GE -20 AND lon LE 45) # (lat GE -90 AND lat LE 40) * ar
; Australia
mask(*,*,5) = (lon GE 90 AND lon LE 180) # (lat GE -90 AND lat LE 0) * ar

out = fltarr(170, 6)
n = 0
.r
for j = 23, 0, -1 do begin
   for i = 0, 35 do begin
      if land[i, j] then begin
         out[n, *] = mask[i, j, *]
         n = n+1
      endif
   endfor
endfor
end

openw, 1, '/scratch/local1/adbethy/idlroutines/regions/regions_areas_TM2grid.dat'
printf, 1, out
close, 1


