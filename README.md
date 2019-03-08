# BETHY-SCOPE-Interactive-Phenology
BETHY land surface model coupled with the SCOPE photosynthesis and fluorescence model. 

Contributors: Alexander Norton (Uni of Melbourne), Peter Rayner (Uni of Melbourne)

Project Description: 
This repository holds the code and working tree development of the coupled model BETHY-SCOPE. 
The BETHY version is Knorr et al. (2010) with interactive phenology. The SCOPE version is 1.53 
van der Tol et al. (2009,doi:10.5194/bg-6-3109-2009) with some slight modifications and bug fixes in the biochemical (chemical.f90) code, coupled similarly to Koffi et al. 
(2015,doi:10.5194/bg-12-4067-2015). 

The project is set up to run the coupled BETHY-SCOPE model to simulate solar-induced fluorescence
(SIF) and compare with satellite SIF observations. We also construct a jacobian matrix from the 
model to determine information content of SIF observations in the context of the model parameters 
and gross primary production (GPP). 
