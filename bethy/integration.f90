!function Fout = integration(F,choice,Ps)
!global nl nlazi lidf nli

SUBROUTINE integration (npt,F,choice,Ps,lidf,Fout)

USE fluo_param, ONLY : nl,nlazi,nli 

!% author: Dr. ir. Christiaan van der Tol (tol@itc.nl)
!% date:     7   December 2007
!% update:   11  February 2008 made modular (Joris Timmermans)
!%
!% function [F2,F3] = F1tot(F,choice,Ps)
!% calculates the layer average and the canopy average of leaf properties
!% per layer, per leaf angle and per leaf azimuth (36)
!%
!% Input:
!%   F       input matrix (3D)   [nli, nlazi,nl]
!%   choice  integration method  'angles'            : integration over leaf angles
!%                               'angles_and_layers' : integration over leaf layers and leaf angles
!%   Ps      fraction sunlit per layer [nl]
!%
!% Output:
!%   Fout    in case of choice = 'angles': [nl]
!%           in case of choice = 'angles_and_layers': [1]


IMPLICIT NONE 

! Input variables 
INTEGER ,INTENT(IN)                              :: npt 
CHARACTER(len=30),INTENT(IN)                     :: choice 
REAL, INTENT(IN), DIMENSION(nli,nlazi,nl)        :: F 
REAL, INTENT(IN), DIMENSION(nl)                  :: Ps 
REAL, INTENT(IN), DIMENSION(nli)                 :: lidf

! Output variables 
REAL, INTENT(OUT), DIMENSION(npt)                :: Fout 

! Local variables 
REAL, DIMENSION(nli,nlazi,nl)                    :: Fout0 
REAL, DIMENSION(nlazi,nl)                        :: Fout1
REAL, DIMENSION(nl)                              :: Fout2
INTEGER                                          :: j 

!Fout = zeros(nli, nlazi,nl);




SELECT CASE (trim(choice))
!%% integration over leaf angles
! Not test for the moment ... EK 5/03/2013
    CASE ('angles')
        print*,'integrating over angles'
        DO j =1, nli 
             Fout0(j,:,:)    = F(j,:,:)*lidf(j)       !  % [nli, nlazi,nl]
        END DO  

        print*,'Fout0:',shape(Fout0)
    !!!    Fout            = permute(Fout,[3 1 2])  !  % [nl]
        Fout1 = sum( Fout0, dim=1 )
        print*,'Fout1:',shape(Fout1)
        Fout2 = sum( Fout1, dim=1 )
        print*,'Fout2:',shape(Fout2)
        DO j=1,nl
            Fout(j) = Fout2(j)/nlazi
        END DO
!        DO j=1,nl 
!             Fout(j) =  sum(Fout0(:,:,j))/nlazi   !  [1,1,nl]
!        END DO    

        print*,'Fout:',shape(Fout)
!%% integration over layers only         
   CASE ('layers')
        !%not implemented
!%% integration over both leaf angles and layers       
   CASE ('angles_and_layers')
        DO j = 1,nli
             Fout0(j,:,:)      = F(j,:,:)*lidf(j)
       !      print*, ' Fout0 ',  Fout0(j,:,:)
        END DO
        
        DO j = 1,nl
             Fout0(:,:,j)    = Fout0(:,:,j)*Ps(j)
       END DO
        Fout                = sum(Fout0)/nlazi/nl
     
END SELECT 

END SUBROUTINE 
