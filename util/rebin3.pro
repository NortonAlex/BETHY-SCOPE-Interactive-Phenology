FUNCTION REBIN3, a, m, dim, missing = missing,  keepsum = keepsum

; rebin array 'a' along dimension number 'dim' to 'm' elements
; equally spaced intervals
; conserves the sum over all entries

s = SIZE (a)
IF s(0) EQ 0 OR s(0) GT 6 THEN $
    Message, 'Array must have 1 to 6 dimensions.'
IF N_ELEMENTS (dim) EQ 0 THEN dim = 1
dim = dim > 1 < s(0)
n = s[dim]

nf = 1l & FOR i = 1, dim - 1 DO nf = nf * s[i]
FOR i = dim + 1, s[0] DO nf = nf * s[i]

ai = s[0]
id = dim-1

if ai gt 2 and dim ne ai then begin
   o = intarr(ai)
   if dim gt 1 then begin
      o[0:id-1] = indgen(id)
      o[id:ai-2] = indgen(ai-id-1)+dim
      o[ai-1] = id
   endif else begin
      o[0] = ai
      o[ai-1] = 0
      o[id:ai-2] = indgen(ai-id-1)+dim
   endelse
   t = FLOAT (FLOAT(transpose(a, o)), 0, nf, n)
endif else begin
   t = FLOAT (FLOAT(a), 0, nf, n)
   if ai eq 2 and dim eq 1 then t = FLOAT (FLOAT(transpose(a)), 0, nf, n)
endelse

b = FLTARR (nf, m)
d = fltarr(nf)

FOR j = 0, m - 1 DO BEGIN
   b[*,j] = 0.
   i0 = FLOAT(j) * n / m
   i1 = FLOAT (j + 1) * n / m
   FOR i = FIX (i0), FIX (i1) < (n - 1) DO begin
      for c = 0l, nf-1 do begin
         if keyword_set(missing) then begin
            if t[c, i] ne missing then begin
               b[c, j] = b[c, j] + (MIN ([i1, i+1]) - MAX ([i0, i])) * t[c, i]
               d[c] = d[c] + (MIN ([i1, i+1]) - MAX ([i0, i]))
            endif 
         endif else begin
            b[c, j] = b[c, j] + (MIN ([i1, i+1]) - MAX ([i0, i])) * t[c, i]
            d[c] = d[c] + (MIN ([i1, i+1]) - MAX ([i0, i]))
         endelse         
      endfor
   endfor
   h = where(d eq 0.)
   If keyword_set(missing) and h[0] ne -1 then  b[h, j] = missing
   if h[0] ne -1 then d[h] = 1.
   if not keyword_set(keepsum) then b[*, j] = b[*, j]/d[*]
   d[*] = 0.

ENDFOR


if ai gt 2 then begin
   r = intarr(ai)
   r[id] = ai-1
   if dim eq 1 then begin
      r[1:ai-1] = indgen(ai-1)
   endif else begin
      r[0:id-1] = indgen(id)
      if dim ne ai then r[id+1:ai-1]=indgen(ai-dim)+id
   endelse
   if dim eq ai then o = indgen(ai)
   o = o + 1
endif

s[dim] = m
CASE ai OF
   1: RETURN, FLOAT (b, 0, s[1])

   2: begin
      if dim eq 1 then RETURN, TRANSPOSE (FLOAT (b, 0, s[2], s[1]))
      RETURN, FLOAT (b, 0, s[1], s[2])
   end

   3: begin
      h = FLOAT (b, 0, s[o[0]], s[o[1]], s[o[2]])
      return,  transpose(h, r)
   end

   4: begin
      h = FLOAT (b, 0, s[o[0]], s[o[1]], s[o[2]], s[o[3]])
      return,  transpose(h, r)
   end

   5: begin
      h = FLOAT (b, 0, s[o[0]], s[o[1]], s[o[2]], s[o[3]], s[o[4]])
      return,  transpose(h, r)
   end

   6: begin
      h = FLOAT (b, 0, s[o[0]], s[o[1]], s[o[2]], s[o[3]], s[o[4]], s[o[5]])
      return,  transpose(h, r)
   end

ENDCASE
RETURN, b
 
END
