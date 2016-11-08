PRO tvf, image, INTERP = interp, BAR = bar, FORMAT = format, NLABELS = n3, ANNOTATION_OFF = noan, $
                REVERSE = rev, MAP = map, COLOR = color, GLINETHICK = glinethick, COMPRESS = comp, $
                AITOFF = aitoff, MOLLWEIDE = mollweide, AZIMUTHAL = azimuthal, STEREOGRAPHIC = stereographic, $
                GNOMIC = gnomic, LAMBERT = lambert, HAMMER=hammer, MERCATOR = mercator, ORTHOGRAPHIC = orthographic, SAT_P = sat_p, $
                ISOTROPIC = iso, HIRES = hires, COUNTRIES = countries, RIVERS = rivers, $
                LATDEL = latdel, LONDEL = londel, TITLE = title, $
                LAT0 = lat0, LON0 = lon0, LIMIT = limit, HORIZON = horizon, $
                LATMIN = latmin, LATMAX = latmax, LONMIN = lonmin, LONMAX = lonmax, $
                GIF_FILE = gif_file,TIFF_FILE = tiff_file, FILE = file, PIXEL = pixel, $
                minscale = minscale, maxscale = maxscale, MISSING = missing, CHARSIZE = chs
                
; as tvscl, but fit into window
; use of color table: ignore last entry for color scale
;                     use last entry as default plot color

; lat0    latitude of centre of map
; lon0    longitude of centre of map
; limit   borders of image in the form [Latmin, Lonmin, Latmax, Lonmax]
; latmin  latitude corresponding to bottom of image
; latmax  latitude corresponding to top of image
; lonmin  longitude corresponding to left edge of image
; lonmax  longitude corresponding to right edge of image

; gif_file, tiff_file  when a file name is specified with this keyword, the resulting
;                      TV image is written in GIF/TIFF format after displaying
; maxscale  highest value for colour table
; minscale  lowest value for colour table
; missing   value for missing data

scalesz = !d.table_size - 2B
top_color = !d.table_size - 1B
scale0 = 1B
IF KEYWORD_SET (file) THEN BEGIN
  SET_PLOT, 'PS'
  DEVICE, /COLOR, /ENCAPS, FILE = file
ENDIF
p0 = [0,0]
IF KEYWORD_SET (bar) THEN p0 = [0, FIX (0.10 * !d.y_size + 0.5)]
topy = 1.
IF KEYWORD_SET (title) THEN topy = .91
xs = !d.x_size - p0(0) & ys = FIX (!d.y_size * topy + 0.5) - p0(1)
;IF NOT KEYWORD_SET (color) THEN color = 255
IF NOT KEYWORD_SET (color) THEN color = top_color
IF NOT KEYWORD_SET (chs) THEN chs = 1.
IF N_ELEMENTS (lat0) EQ 0 THEN lat0 = 0.
IF N_ELEMENTS (lon0) EQ 0 THEN lon0 = 0.
s = SIZE (image)
; number of leading dimensions equal to one:
i1 = WHERE (s(1:*) GT 1) & i1 = i1(0)
IF N_ELEMENTS (latmin) EQ 0 THEN latmin =  -90 + 180. / s(2+i1) / 2.
IF N_ELEMENTS (latmax) EQ 0 THEN latmax =   90 - 180. / s(2+i1) / 2.
IF N_ELEMENTS (lonmin) EQ 0 THEN lonmin = -180 + 360. / s(1+i1) / 2.
IF N_ELEMENTS (lonmax) EQ 0 THEN lonmax = lonmin + 360. - 360. / s(1+i1)
IF KEYWORD_SET (missing) THEN BEGIN
  ; produce list of valid points
  l = WHERE (image NE missing)
  IF l(0) NE -1 THEN lo = min (image(l)) ELSE lo = min (image)
  IF l(0) NE -1 THEN hi = max (image(l)) ELSE hi = max (image)
ENDIF ELSE BEGIN
  lo = min (image)
  hi = max (image)
ENDELSE
IF N_ELEMENTS (minscale) EQ 0 THEN minscale = lo
IF N_ELEMENTS (maxscale) EQ 0 THEN maxscale = hi

;ERASE, 0B

IF s(0) GT 2 - i1 AND s(1) NE 1 THEN BEGIN
  IF NOT KEYWORD_SET (comp) THEN comp = 4
  n = s(s(0)+2) / s(1+i1) / s(2+i1)
  ny = (n-1) / 3 + 1
  nx = 3
  a = FLOAT (image, 0, s(1+i1), s(2+i1), n)
  IF KEYWORD_SET (rev) THEN a = REVERSE (a, 2)
  IF NOT KEYWORD_SET (map) THEN BEGIN
    a = CONGRID (a, xs/nx, ys/ny, n, INTERP = interp)
  ENDIF ELSE BEGIN
    IF KEYWORD_SET (missing) THEN l = WHERE (a EQ missing) ELSE l = -1
    a = BYTSCL (a, MIN = minscale, MAX = maxscale, top = scalesz - 1) + scale0
    IF l(0) NE -1 THEN a(l) = 0B
  ENDELSE
  FOR j = 0, ny - 1 DO FOR i = 0, nx - 1 DO BEGIN
    k = i+j*3 & x = p0(0)+i*FLOAT(xs)/nx & y = p0(1)+(ny-1-j)*FLOAT(ys)/ny
    IF k LT n THEN BEGIN
      IF NOT KEYWORD_SET (map) THEN BEGIN
        ; rescale image in bytes, 0 and highest color index are left out
        img = BYTSCL (a(*,*,k), MIN = minscale, MAX = maxscale, top = scalesz - 1) + scale0
        ; treat missing data
        IF KEYWORD_SET (missing) THEN l = WHERE (a(*,*,k) EQ missing) ELSE l = -1
        IF l(0) NE -1 THEN img(l) = 0B
        ERASE, 0B
        TV, img, x, y
        IF KEYWORD_SET (title) THEN xyouts, 0.5, topy+.05, title, ALIGNMENT = 0.5, CHARSIZE = chs*1.2, CHARTHICK=chs, /NORMAL, COLOR=1B
      ENDIF ELSE BEGIN
        MAP_SET, /NOBORDER, POSITION = [x/!d.x_size, y/!d.y_size, $
          (x+xs/nx)/!d.x_size, (y+ys/ny)/!d.y_size], /ADVANCE, $
          AITOFF = aitoff, MOLLWEIDE = mollweide, AZIMUTHAL = azimuthal, STEREOGRAPHIC = stereographic, $
          GNOMIC = gnomic, LAMBERT = lambert, HAMMER = hammer, MERCATOR = mercator, ORTHOGRAPHIC = orthographic, SAT_P = sat_p, $
          lat0, lon0, LIMIT = limit, HORIZON = horizon, ISOTROPIC = iso, $
          CHARSIZE = chs*1.3, GLINETHICK = glinethick
        img = MAP_IMAGE (a(*,*,k), xpos, ypos, xsz, ysz, $
              latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, missing=missing)
        ERASE, 0B
        TV, img, xpos, ypos, XSIZE = xsz, YSIZE = ysz
        MAP_CONTINENTS, HIRES=hires, /COASTS, COUNTRIES=countries, RIVERS=rivers, COLOR = color, MLINETHICK = glinethick
        MAP_GRID, COLOR = color, LATDEL = latdel, LONDEL = londel, GLINETHICK = glinethick
        IF KEYWORD_SET (title) THEN xyouts, 0.5, topy+.05, title, ALIGNMENT = 0.5, CHARSIZE = chs*1.2, CHARTHICK=chs, /NORMAL, COLOR=1B
      ENDELSE
      XYOUTS, x, y, STRCOMPRESS (k, /REM), ALIGN = 0., /DEVICE, COLOR = color, CHARSIZE = chs, CHARTHICK = chs
    ENDIF
  ENDFOR
ENDIF ELSE BEGIN
  IF NOT KEYWORD_SET (map) THEN BEGIN
    x = p0(0) & y = p0(1) & xsz = xs & ysz = ys
    IF i1 EQ 0 THEN $
      a = CONGRID (image, xsz, ysz, INTERP = interp) $
    ELSE IF i1 EQ 1 THEN $
      a = CONGRID (image, 1, xsz, ysz, INTERP = interp) $
    ELSE MESSAGE, 'TVF: Input array must have 2 or 3 dimensions'
    IF KEYWORD_SET (rev) THEN a = REVERSE (a, 2+i1)
    ; rescale image in bytes, 0 and highest color index are left out
    img = BYTSCL (a, MIN = minscale, MAX = maxscale, top = scalesz - 1) + scale0
    ; treat missing data
    IF KEYWORD_SET (missing) THEN l = WHERE (a EQ missing) ELSE l = -1
    IF l(0) NE -1 THEN img(l) = 0B
    ERASE, 0B
    TV, img, x, y
    IF KEYWORD_SET (title) THEN xyouts, 0.5, topy-.005, title, ALIGNMENT = 0.5, CHARSIZE = chs*1.3, CHARTHICK = chs, /NORMAL, COLOR=1B
  ENDIF ELSE BEGIN
    IF NOT KEYWORD_SET (comp) THEN comp = 1
    MAP_SET, /NOBORDER, POSITION = [FLOAT(p0(0))/xs, FLOAT(p0(1))/ys, 1., topy], $
       AITOFF = aitoff, MOLLWEIDE = mollweide, AZIMUTHAL = azimuthal, STEREOGRAPHIC = stereographic, $
       GNOMIC = gnomic, LAMBERT = lambert, HAMMER = hammer, MERCATOR = mercator, ORTHOGRAPHIC = orthographic, $
       ISOTROPIC = iso, GLINETHICK = glinethick, $
       lat0, lon0, LIMIT = limit, HORIZON = horizon
    a = BYTSCL (image, MIN = minscale, MAX = maxscale, top = scalesz - 1) + scale0
    ; treat missing data
    IF KEYWORD_SET (missing) THEN l = WHERE (image EQ missing) ELSE l = -1
    IF l(0) NE -1 THEN a(l) = 0B
    img = MAP_IMAGE (a, x, y, xsz, ysz, COMPRESS = comp, $
          latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax)
    ERASE, 0B
    TV, img, x, y, XSIZE = xsz, YSIZE = ysz
    MAP_CONTINENTS, HIRES=hires, /COASTS, COUNTRIES=countries, RIVERS=rivers, COLOR = color, MLINETHICK = glinethick
    MAP_GRID, COLOR = color, LATDEL = latdel, LONDEL = londel, GLINETHICK = glinethick
    IF KEYWORD_SET (title) THEN $
      xyouts, 0.5, topy+.05, title, ALIGNMENT = 0.5, CHARSIZE = chs*1.2, CHARTHICK=chs, /NORMAL, COLOR=1B
  ENDELSE
ENDELSE

IF KEYWORD_SET (bar) THEN BEGIN
  p1 =  FIX (CONVERT_COORD ([0.35,0.035], /NORMAL, /TO_DEVICE) + 0.5) ; bar, lower left
  p2 =  FIX (CONVERT_COORD ([0.65,0.060], /NORMAL, /TO_DEVICE) + 0.5) ; bar, upper right
  p3 =  FIX (CONVERT_COORD ([0.35,0.005], /NORMAL, /TO_DEVICE) + 0.5) ; annotation, lower left
  IF NOT KEYWORD_SET (n3) THEN n3 = 3
  a = BYTSCL (REBIN (FINDGEN (p2(0)-p1(0)), p2(0)-p1(0), p2(1)-p1(1)), top = scalesz - 1) + scale0
  TV, a, p1(0), p1(1)
  PLOTS, [p1(0), p2(0), p2(0), p1(0), p1(0)], [p1(1), p1(1), p2(1), p2(1), p1(1)], /DEVICE, COLOR = color
  hi = maxscale & lo = minscale
  lbl = lo + INDGEN (n3) * FLOAT (hi - lo) / (n3 - 1)
  IF NOT KEYWORD_SET (format) THEN format = '(E9.2)'
  IF KEYWORD_SET (format) THEN lbl = STRING (lbl, FORMAT = format)
  IF NOT KEYWORD_SET (noan) THEN $
  FOR i = 0, n3 - 1 DO XYOUTS, p1(0) + i*(p2(0)-p1(0))/(n3-1), p3(1), CHARSIZE = chs, CHARTHICK = chs, $
    STRCOMPRESS (lbl(i), /REM), ALIGN = 0.5, /DEVICE, COLOR = color
  FOR i = 0, n3 - 1 DO PLOTS, [p1(0) + i*(p2(0)-p1(0))/(n3-1), p1(0) + i*(p2(0)-p1(0))/(n3-1)], $
    [0.2*p3(1)+0.8*p1(1),p1(1)], /DEVICE, COLOR = color
ENDIF
  
IF KEYWORD_SET (pixel) THEN BEGIN
  hi = maxscale & lo = minscale
  img = TVRD ()
  RDPIX, FLOAT(img) / MAX (img) * (hi - lo) + lo
ENDIF

IF KEYWORD_SET (file) THEN BEGIN
  DEVICE, /CLOSE
  SET_PLOT, 'X'
ENDIF

IF KEYWORD_SET (gif_file) THEN BEGIN
  TVLCT, r, g, b, /GET
  PRINT, '* writing to GIF file: ' + gif_file + ' ...'
  WRITE_GIF, gif_file, TVRD (), r, g, b
ENDIF

IF KEYWORD_SET (tiff_file) THEN BEGIN
  PRINT, '* writing to TIFF file: ' + tiff_file + ' ...'
  DEVICE, GET_DECOMPOSED = currentDecomposed
  IF currentDecomposed EQ 1 THEN BEGIN
    TVLCT, r, g, b, /GET    WRITE_TIFF, tiff_file, reverse(TVRD (),2), $
      red=r, green=g, blue=b, compression=2  ENDIF ELSE BEGIN
    store_image =  COLOR_QUAN (reverse(TVRD (true=1),3), 1, r, g, b)
    WRITE_TIFF, tiff_file, store_image, $
      red=r, green=g, blue=b, compression=2  ENDELSE
ENDIF

END
