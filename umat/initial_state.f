	  SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      implicit none
C
      REAL(8) STATEV(NSTATV),COORDS(NCRDS)            
      integer NOEL,NCRDS,NPT,LAYER,KSPT,NSTATV
	  
      statev(:)   = 0.0d0
      statev(1)   = 0.685   !epor -make it similar to epor at nodes for cpe4p elements
      statev(2)   = 0.0   ! pwater
      

      RETURN
      END


 

   
