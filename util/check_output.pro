PRO check_output, in, lat, lon, m, x, y, field = field, format = format, overmax = overmax, undermin = undermin, mean=mean

; USAGE:
; from sandbox do
; cd ./util
; idl
; check_output
; OR
; checkoutput, field = 'fpar'
; ETC.
 
IF NOT KEYWORD_SET (field) THEN field = 'lai'
field = STRCOMPRESS (field,/REM)

path = '../output/'
device,decompose=0
loadct,40
tvlct, 255B, 255B, 255B
miss = -1e6
loadm,path+field+'.dat',gr=2.,/center,year=-1,in,oc=miss,format=format

s = size (in) & nx = s[1] & ny = s[2] & nm = s[3]
lon_list = (findgen (nx) + 0.5) / nx * 360. - 180.
lat_list = (findgen (ny) + 0.5) / ny * 180. - 90.
t = (findgen(nm) + 0.5) / 12.
nyr = nm / 12

tvf,total(in,3)/nm,missing=miss,/map,/bar,tif=field+'.tif',title=STRCOMPRESS(nyr,/REM)+'-year mean '+field+' from CCDAS-BETHY'

IF KEYWORD_SET (mean) THEN BEGIN
  in0 = total(in,3)/nm
  IF N_ELEMENTS (overmax) NE 0 THEN l = where (in0 GT overmax AND in0 NE miss)
  IF N_ELEMENTS (undermin) NE 0 THEN l = where (in0 LT undermin AND in0 NE miss)
ENDIF ELSE BEGIN
  IF N_ELEMENTS (overmax) NE 0 THEN l = where (in GT overmax AND in NE miss)
  IF N_ELEMENTS (undermin) NE 0 THEN l = where (in LT undermin AND in NE miss)
ENDELSE

IF N_ELEMENTS (overmax) NE 0 OR N_ELEMENTS (undermin) NE 0 THEN BEGIN
  m = l/(1L*nx*ny)
  y = (l-m*(nx*ny)) / nx
  x = l - m*(nx*ny) - y*nx
  lon = lon_list[x]
  lat = lat_list[y]
  STOP
ENDIF

;tvf,total(in,3)/nm<0,missing=miss,/map,/bar

END
                                                                                             
