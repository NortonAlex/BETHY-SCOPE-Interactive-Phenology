;***************************************************************************
;***************************************************************************
; IDL functions:
;
;    LoadProfileData_RedRes() 
;                To load FAPAR Product containing time
;                composite, time profile data.
;    ScaleProfileData()
;                To correctly scale the data loaded with
;                LoadProfileData_RedRes().
;
; Internal functions (for program use only):
;    _ConvertDate()
;                To convert a file embedded date (year-month-day) to [day of
;                year, year] format.
;    _FileStruct()
;                Depending on keywords, either: 
;                   parse the file to extract information to
;                   establish the data structure
;                or: 
;                   load the data structure from the file.
;
;***************************************************************************

;---------- Internal Functions -----------
function _ConvertDate, dateString
   d = strsplit(dateString, '-', /extract)
   year = fix(d[0])
   month = fix(d[1])
   day = fix(d[2])
   day = julday(month, day, year) - julday(1, 1, year) + 1
   return, {day:day, year:year}
end

function _FileStruct,  fName, dataStruct=dataStruct, MERIS=meris, TC_MERIS=tc_meris

		varsFlag = intarr(16) 
		varsType = [1,1,1,1,1,1,1,2,2,2,2,4,4,4,4,1*3]
nCols = 0l
nRows = 0l
fieldStr = ''
fromStr = ''
toStr = ''
varName = ''
dateCount = -1

openr, unit, fName, /get_lun
while( not(EOF(unit))) do begin
	found = -1
	readf, unit, fieldStr
	case fieldStr of
	'Columns'  : begin
		readf, unit, nCols
		dateCount = dateCount + 1
	end
	'Rows'     : readf, unit, nRows
	'Date From': begin
		readf, unit, fromStr
		if(keyword_set(dataStruct)) then begin
			dT = _ConvertDate(fromStr)
			dataStruct[dateCount].startDay = dT.day
			dataStruct[dateCount].year = dT.year
		endif
	end
	'Date To'  : begin
		readf, unit, toStr
		if(keyword_set(dataStruct)) then begin
			dT = _ConvertDate(toStr)
			dataStruct[dateCount].endDay = dT.day
		endif
	end
	'Variable' : begin
		readf, unit, varName
			case varName of
			'Nearest Neighbour:flag'        :  found = 0
			'flag'        :  found = 0
			'Nearest Neighbour:MGVI'       :  found = 1
			'MGVI'       :  found = 1
			'Nearest Neighbour:dMGVI'      :  found = 2
			'dMGVI'      :  found = 2
			'Nearest Neighbour:sd_MGVI'    :  found = 3
			'sd_MGVI'    :  found = 3
			'Nearest Neighbour:nb_MGVI'    :  found = 4
			'nb_MGVI'    :  found = 4
			'Nearest Neighbour:BRF_Rec_Red'   :  found = 5
			'BRF_Rec_Red'   :  found = 5
			'Nearest Neighbour:BRF_Rec_Nir'   :  found = 6
			'BRF_Rec_Nir'   :  found = 6
			'Nearest Neighbour:norm_surf_reflec_2' :  found = 7
			'norm_surf_reflec_2' :  found = 7
			'Nearest Neighbour:norm_surf_reflec_5' :  found = 8
			'norm_surf_reflec_5' :  found = 8
			'Nearest Neighbour:norm_surf_reflec_8' :  found = 9
			'norm_surf_reflec_8' :  found = 9
			'Nearest Neighbour:norm_surf_reflec_13' :  found = 10
			'norm_surf_reflec_13' :  found = 10
			'Nearest Neighbour:solar_zenith' :  found = 11
			'solar_zenith' :  found = 11
			'Nearest Neighbour:view_zenith'  :  found = 12
			'view_zenith'  :  found = 12
			'Nearest Neighbour:solar_azimuth'  :  found = 13
			'solar_azimuth'  :  found = 13
			'Nearest Neighbour:view_azimuth' :  found = 14
			'view_azimuth' :  found = 14
			'Nearest Neighbour:Nearest Neighbour:flag_ass_pix.flags'   :  found = 15
			'Nearest Neighbour:flag_ass_pix.flags'   :  found = 15
			else: print,'bizarre TC_MERIS '
			;+varName
			endcase
		if(found ge 0) then begin
			if( keyword_set(dataStruct)) then begin
				tData = dataStruct[dateCount].(found)
				readu, unit, tData,TRANSFER_COUNT=tc
;				if(size(tData, /type) eq 2) then byteorder, tData, /SWAP_IF_LITTLE_ENDIAN ;wrong
				dataStruct[dateCount].(found) =  tData
			endif else begin
				len = nRows * nCols
				inc = varsType(found) * len
				varsFlag[found] = varsFlag[found] + 1

				point_lun, -unit, currPos 
				currPos = currPos + inc
				point_lun, unit, currPos
				
				;print, 'Variable name = '+varname 
			endelse
		endif
	; end case 'Variable'
	end
	else: print,'bizarre case '
	endcase
    endwhile

    IF (keyword_set(dataStruct)) THEN BEGIN
	    return, 1 
    ENDIF ELSE BEGIN
	    xL = (varsFlag gt 0) * nCols
	    yL = (varsFlag gt 0) * nRows
	    xL = xL > 1
	    yL = yL > 1
	    
	    dStruct = {day:0b, month:0b, year:0, yearDay:0}
		    tStruct = {$
			    flag:bytarr(xL[0], yL[0]), $
			    fapar:bytarr(xL[1], yL[1]), $
			    daySelected:bytarr(xL[2],yL[2]), $
			    sd_fapar:bytarr(xL[3],yL[3]), $
			    n_fapar:bytarr(xL[4],yL[4]), $
			    brf_rec_r:bytarr(xL[5],yL[5]),$
			    brf_rec_n:bytarr(xL[6],yL[6]), $
			    brf_norm_2:uintarr(xL[7],yL[7]), $
			    brf_norm_8:uintarr(xL[8],yL[8]),$
			    brf_norm_5:uintarr(xL[9],yL[9]), $
			    brf_norm_13:uintarr(xL[10],yL[10]), $
			    sun_zenith_MERIS:lonarr(xL[11],yL[11]), $
			    sat_zenith_MERIS:lonarr(xL[12],yL[12]),$
			    sun_azimuth_MERIS:lonarr(xL[13],yL[13]), $
			    sat_azimuth_MERIS:lonarr(xL[14],yL[14]),$
			    flag_MERIS:bytarr(xL[15],yL[15],3), $
			    startDay:0, endDay:0, year:0}
	    nStructs = max(varsFlag)
	    close, unit
	    free_lun, unit
	    return, replicate(tStruct, nStructs)
    ENDELSE
END
    
function _FloatStruct, intStruct
    nTags = n_tags(intStruct)
    dL = lonarr(2,nTags)
    for i=0, nTAgs - 1 do dL[*,i] = size(intStruct[0].(i), /dimensions)
    dL = dL > 1
    return, {$
	    flag:bytarr(dL[0,0], dL[1,0]),$ 
    fapar:fltarr(dL[0,1], dL[1,1]), $
	    daySelected:bytarr(dL[0,2],dL[1,2]),$
	    sd_fapar:fltarr(dL[0,3],dL[1,3]), n_fapar:bytarr(dL[0,4],dL[1,4]),$ 
    brf_rec_r:fltarr(dL[0,5],dL[1,5]),$
	    brf_rec_n:fltarr(dL[0,6],dL[1,6]), brf_norm_2:fltarr(dL[0,7],dL[1,7]), brf_norm_8:fltarr(dL[0,8],dL[1,8]),$
	    brf_norm_5:fltarr(dL[0,9],dL[1,9]), brf_norm_13:fltarr(dL[0,10],dL[1,10]), $
	    sun_zenith:fltarr(dL[0,11],dL[1,11]), sat_zenith:fltarr(dL[0,12],dL[1,12]), $
	    sun_azimuth:fltarr(dL[0,13],dL[1,13]), sat_azimuth:fltarr(dL[0,14],dL[1,14]),$ 
    startDay:0, endDay:0, year:0}
    end
    
    ;-------------- End internal functions. --------------------

;***************************************************************************
;Function to load a JRC-FAPAR time profile,
;data block.
;
; Usage:
;   siteProfile = LoadProfileData_RedRes( file )

; Where:
;   file is the full path and name of the raw data file to be loaded.
;   
; Keywords:
;   None
;
; Returns:
;   A vector of structures, where the vector dimension corresponds to
;   time, and the structure corresponds to a single, time composite
;   product for a particular time.
;   It is possible that not all fields were supplied in a file, under
;   these circumstances the array size for these particular fields will be 1.
;
; Structure:
;   The basic struture element is:
;   {
;       flag:bytarr(xLen, ylen), 
;       fapar:bytarr(xLen, ylen), 
;       daySelected:bytarr(xLen,yLen), 
;       sd_fapar:bytarr(xLen,yLen), 
;       n_fapar:bytarr(xLen,yLen)
;       brf_rec_r:intarr(xLen,yLen), 
;       brf_rec_n:intarr(xLen,yLen), 
;       brf_TOA_443:intarr(xLen,yLen), 
;       brf_TOA_670:intarr(xLen,yLen), 
;       brf_TOA_865:intarr(xLen,yLen), 
;       sun_zenith:intarr(xLen,yLen), 
;       sat_zenith:intarr(xLen,yLen), 
;       rel_azimuth:intarr(xLen,yLen), 
;       startDay:0, 
;       endDay:0, 
;       year:0, 
;   }
;
; Implementation
;
;   Examine the input file structure and create an appropriate IDL data
;   structure.
;   Not all fields need be in the file.
;   Load the data from the file.
;   On input, check if the data need to be byte swapped, and do so if
;   necessary.
;   

;
function LoadProfileData_RedRes,  fName 
   dataStruct = _FileStruct(fName)
   status     = _FileStruct(fName, dataStruct=dataStruct)
  return, dataStruct
end
;***************************************************************************
;***************************************************************************
