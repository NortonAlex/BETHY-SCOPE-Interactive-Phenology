FUNCTION bethy2tm, a,  miss = miss,  keep = keep

; CONVERT TO TM2 COARSE-GRID FORMAT
; INPUT UNITS: no units, 180x90 array
; OUTPUT UNITS: no units, 36x24 array

s = SIZE (a)
nx = s[1] & ny = s[2] 
ni = 36 & nj = 24
if s[0] eq 3 then nk = s[3]
if s[0] eq 4 then begin
   nk = s[3]
   nl = s[4]
endif

;bring to TM2 resolution in kgC / grid cell and year
hj = 2*nj-2
n2 = a
if keyword_set(miss) then begin
   if keyword_set(keep) then begin
      n2 = rebin3(n2, ni, 1, missing = miss, keepsum = keep)
      n2 = rebin3(n2, 2*nj-2, 2, missing = miss, keepsum = keep) 
   endif else begin
      n2 = rebin3(n2, ni, 1, missing = miss)
      n2 = rebin3(n2, 2*nj-2, 2, missing = miss) 
   endelse
endif else begin
   if keyword_set(keep) then begin
      n2 = rebin3(n2, ni, 1, keepsum = keep)
      n2 = rebin3(n2, 2*nj-2, 2, keepsum = keep) 
   endif else begin
      n2 = rebin3(n2, ni, 1)
      n2 = rebin3(n2, 2*nj-2, 2) 
   endelse
endelse

CASE s[0] OF
   2: begin
         net0 = FLTARR (ni, nj)
         ; polar grid value stored in first element only (O.K.)
         net0[0, 0] = n2[0, 0]  ; SOUTH POLE
         net0[0, nj-1] = n2[0, 2*nj-3] ; NORTH POLE
         if keyword_set(keep) then begin
            net0[*, 1:nj-2] = REBIN3 (n2[*, 1:2*nj-4], nj-2, 2, keepsum = keep) 
         endif else begin
            net0[*, 1:nj-2] = REBIN (n2[*, 1:2*nj-4], ni, nj-2) 
         endelse
      end

   3: begin
         net0 = FLTARR (ni, nj, nk)
         ; polar grid value stored in first element only (O.K.)
         net0[0, 0, *] = n2[0, 0, *] ; SOUTH POLE
         net0[0, nj-1, *] = n2[0, 2*nj-3, *] ; NORTH POLE
         if keyword_set(keep) then begin
            net0[*, 1:nj-2, *] = REBIN3 (n2[*, 1:2*nj-4, *], nj-2, 2, keepsum = keep) 
         endif else begin
            net0[*, 1:nj-2, *] = REBIN (n2[*, 1:2*nj-4, *], ni, nj-2, nk) 
         endelse
      end

   4: begin
         net0 = FLTARR (ni, nj, nk, nl)
         ; polar grid value stored in first element only (O.K.)
         net0[0, 0, *, *] = n2[0, 0, *, *] ; SOUTH POLE
         net0[0, nj-1, *, *] = n2[0, 2*nj-3, *, *] ; NORTH POLE
         if keyword_set(keep) then begin
            net0[*, 1:nj-2, *, *] = REBIN3 (n2[*, 1:2*nj-4, *, *], nj-2, 2, keepsum = keep) 
         endif else begin
            net0[*, 1:nj-2, *, *] = REBIN (n2[*, 1:2*nj-4, *, *], ni, nj-2, nk, nl) 
         endelse
      end
ENDCASE

RETURN, net0

END
