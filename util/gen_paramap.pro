PRO gen_paramap, lores = lores, site = site

; GENERATE PARAMETER MAPPING FILES FOR HIRES, LORES OR LIST OF SITES
; Output:
; File named hires_XXparams.map in output directory ../control_bethy.
; where XX is the number of parameters at a single sub-gridcell.
; Options:
; /lores       generate file lores_XXparams.map for the lores grid, instead of the default
; site=<name>  specify site specific file <name> in same directory as output,
;              e.g. 'site_specs.dat', instead of default; writes to file <name>.map

; SPECIFICATIONS
res = 'hires'
IF KEYWORD_SET (lores) THEN res = 'lores'
;np = 27   ; number of parameters at single model sub-pixel
;np = 31   ; number of parameters at single model sub-pixel
;np = 27   ; number of parameters at single model sub-pixel
np = 28   ; number of parameters at single model sub-pixel
npft = 13 ; maximum number of PFTs
odir = '../control_bethy/'  ; output directory
;
IF NOT KEYWORD_SET (site) THEN BEGIN
  loadl,'gridinfo_'+res+'.txt',d,6,100000
  ng = n_elements (d[0,*])
  pft = fix(reform(d[2,*]))
  openw,u,/get_lun,odir+res+'_'+strcompress(np,/rem)+'params.map'
ENDIF ELSE BEGIN
  ng = 0
  ngmax = 10000
  maxpft_site = 3
  pft = intarr(ngmax)
  pftt = intarr(maxpft_site)
  openr,v,/get_lun,odir+site
  line = ' '
  readf,v,line
  WHILE NOT EOF(v) AND ng LE ngmax-maxpft_site DO BEGIN
    readf,v,lat,lon,elev,pftt,format='(2f7.2,i6,3i3)'
    print,lat,lon,elev,pftt
;    readf,v,pftt,format='(20X,'+strcompress(maxpft_site,/rem)+'I3)'
    n = fix(total(pftt NE 0)) ; ignore sub-gridcells with PFT=0
    pft[ng:ng+n-1] = pftt(where (pftt NE 0))
    ng = ng + n
  ENDWHILE
  pft = pft[0:ng-1]
  free_lun, v
  openw,u,/get_lun,odir+site+'.map'
ENDELSE
;
pftpars = [1, 2, 20] ; PFT depdendent ones among parameters 1...np
offset = intarr (np,npft)
;offset[21,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2] ; parameter 22
;offset[22,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 23
;offset[26,*] = [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1] ; parameter 27
;offset[22,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2] ; parameter 23
;offset[23,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 24
;offset[27,*] = [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1] ; parameter 28
;
;offset[22-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2] ; parameter 22
;offset[23-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 23
;offset[24-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2] ; parameter 24
;offset[25-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 25
;offset[27-1,*] = [1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1] ; parameter 27
;offset[28-1,*] = [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1] ; parameter 28
;                 1  2  3  4  5  6  7  8  9 10 11 12 13
;offset[22-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2] ; parameter 22
;offset[23-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 23
;offset[26-1,*] = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0] ; parameter 26
;offset[27-1,*] = [0, 1, 0, 0, 0, 0, 0, 1, 2, 2, 0, 2, 2] ; parameter 27
;                 1  2  3  4  5  6  7  8  9 10 11 12 13
;offset[22-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2] ; parameter 22
;offset[23-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 23
;offset[27-1,*] = [0, 2, 0, 2, 1, 2, 0, 2, 2, 2, 1, 2, 2] ; parameter 27
;offset[28-1,*] = [0, 1, 0, 0, 0, 0, 0, 1, 2, 2, 0, 2, 2] ; parameter 28
;                 1  2  3  4  5  6  7  8  9 10 11 12 13
;offset[22-1,*] = [0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 3] ; parameter 22
;offset[23-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 23
;offset[27-1,*] = [0, 2, 0, 2, 1, 2, 0, 2, 2, 2, 1, 2, 2] ; parameter 27
;offset[28-1,*] = [0, 1, 0, 0, 0, 0, 0, 1, 2, 2, 0, 2, 2] ; parameter 28
;                 1  2  3  4  5  6  7  8  9 10 11 12 13
;offset[22-1,*] = [0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 3] ; parameter 22
;offset[23-1,*] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0] ; parameter 23
;offset[27-1,*] = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0] ; parameter 27
;offset[28-1,*] = [0, 1, 0, 0, 0, 0, 0, 0, 2, 2, 0, 2, 2] ; parameter 28
;                 1  2  3  4  5  6  7  8  9 10 11 12 13
offset[22-1,*] = [0, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 3, 4] ; parameter 22
offset[23-1,*] = [0, 0, 0, 0, 1, 1, 0, 0, 2, 2, 1, 2, 0] ; parameter 23
offset[27-1,*] = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0] ; parameter 27
offset[28-1,*] = [0, 1, 0, 0, 0, 0, 0, 0, 2, 2, 0, 2, 2] ; parameter 28
FOR i = 0, n_elements(pftpars) - 1 DO offset[pftpars[i]-1,*] = indgen(npft)
;
printf, u, np, ng, format='(I6,I10)'
p = intarr (np)
FOR g=0,ng-1 DO BEGIN
   k = 1
   FOR i=0,np-1 DO BEGIN
      p[i] = k + offset[i,pft[g]-1]
      k = k + 1 + max(offset[i,*])
   ENDFOR
   printf, u, p, format = '(999I6)'
ENDFOR
free_lun, u

END
