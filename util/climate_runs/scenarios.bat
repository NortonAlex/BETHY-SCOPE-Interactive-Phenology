; batch job for generating climate scenarios
;
; compile the code
.rnew ipcc_HadCM3_A1           ; converts HadCM3 scenarios to 2x2 grid
.rnew ipcc_HadCM3_climatology  ; same, but generates climatology
.rnew ipcc_ECHAM5_A1           ; converts ECHAM5 scenarios to 2x2 grid
.rnew ipcc_ECHAM5_climatology  ; same, but generates climatology
.rnew daily_to_monthly         ; generates monthly observations
.rnew daily_climate_parameters ; generates parameters for weather generator
.rnew weather_generator        ; the actual weather generator
;
; sequence for generating weather data for sites, using HadCM
;ipcc_HadCM3_A1
;ipcc_HadCM3_climatology
;daily_to_monthly
;;daily_climate_parameters, /sites
weather_generator, /sites
;
; sequence for generating weather data for sites, using ECHAM5
;ipcc_ECHAM5_A1
;ipcc_ECHAM5_climatology
;daily_to_monthly
;daily_climate_parameters, /sites
weather_generator, /sites, gcm='ECHAM5'
;
; sequence for generating weather data globally, using HadCM
;ipcc_HadCM3_A1
;ipcc_HadCM3_climatology
;daily_to_monthly
;;daily_climate_parameters
weather_generator
;
; sequence for generating weather data globally, using ECHAM5
;ipcc_HadCM3_A1
;ipcc_HadCM3_climatology
;daily_to_monthly
;daily_climate_parameters
weather_generator, gcm='ECHAM5'
;
; sequence for generating weather data globally in lores, using HadCM3
;ipcc_HadCM3_A1
;ipcc_HadCM3_climatology
;daily_to_monthly
;;daily_climate_parameters, /lores
weather_generator, /lores
;
; sequence for generating weather data globally in lores, using ECHAM5
;ipcc_ECHAM5_A1
;ipcc_ECHAM5_climatology
;daily_to_monthly
;daily_climate_parameters, /lores
weather_generator, /lores, gcm='ECHAM5'
;
;
weather_generator
weather_generator, gcm='ECHAM5'
weather_generator, /lores
weather_generator, /lores, gcm='ECHAM5'
weather_generator, /sites
weather_generator, /sites, gcm='ECHAM5'
