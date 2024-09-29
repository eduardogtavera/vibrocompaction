      !---------------------------------------------------------------
      subroutine get_source_position(coord_source,is_active,kstep,kinc,time)
      ! returns the position of the source at a given time.
      ! if is_active=.false. there is no vibration at the source
      use globalVariables
      implicit none
      integer, intent(in) :: kstep
      integer, intent(in) :: kinc 
      real(8), intent(in), dimension(2) ::  time
      real(8), intent(out), dimension(3) :: coord_source
      logical, intent(out) :: is_active
      ! local variables
      !integer, parameter :: activation_step = 2
      real(8), parameter :: t_start = 0 ! velocity of penetration 0.05 m/s
      real(8), parameter :: t_end = 1250 ! penetration + vibration: upwars velocity during densification 0.0066 m/s (150 s per 1 m)
      real(8), parameter :: t_end_penetration = 200.0d0 ! velocity of penetration 0.05 m/s (downwards)
      ! vibration phase
      real(8), parameter :: v_mean_up =  1.0d0/150.0d0  ! mean upwards velocity 0.0066 m/s 
      real(8), parameter :: u_amp =  0.5d0  ! amplitude of oscillation during vibration
      real(8), parameter :: period =  40.0d0  ! period of an oscillation during vibration
      !
      real(8), dimension(3) :: pos_start = [0.0d0, 0.0d0,0.0d0]
      real(8), dimension(3) :: pos_end   = [0.0d0, 0.0d0, 0.0d0]
      real(8) :: alpha
      real(8) :: t_vibration 
	  real(8) :: xtop, ytop
	  integer :: activation_step, kpile						
      
      
      is_active = .False. ! per default, there is no vibration at the source
      coord_source(:) = 0.0d0
	  
	  
	  kpile = inst_sequence(kstep)
	  if (kpile .ne. 0) activation_step = kstep
	  
	  
	  
	  xtop = bhole_top(kpile,1)
	  ytop = bhole_top(kpile,2)
	  pos_start(1:2) = [xtop, ytop]
	  pos_end(1:3) = [xtop, ytop, -9.0d0]
	  
	  

      ! simple example. a source penetrating at a constant vel. in the ground
      if (activation_step .ne. kstep) then 
         return 
      end if
      if ((time(1) .lt. t_start) .or. (t_end .lt. time(1))) then 
         return 
      end if
	  
	  
	  
      is_active = .true.
      if(time(1) .lt. t_end_penetration ) then
         alpha = (time(1)-t_start)/(t_end_penetration-t_start)
         coord_source = pos_start + (pos_end - pos_start)*alpha ! current position of the source
      end if 
	  
      if (time(1) .ge. t_end_penetration) then 
         t_vibration  = time(1) - t_end_penetration		 
		 
         call get_source_position_vibration(coord_source, pos_end, 
     &                              v_mean_up, u_amp, period, t_vibration )
      end if  
	  
      end subroutine get_source_position      
      !---------------------------------------------------------------


      !---------------------------------------------------------------
      subroutine get_source_position_vibration(pos, pos_start, 
     &                         v_mean_up, u_amp, period, curr_t)
      implicit none
      real(8), intent(in), dimension(3) :: pos_start 
      real(8), intent(out), dimension(3) :: pos 
      real(8), intent(in) :: v_mean_up
      real(8), intent(in) :: u_amp
      real(8), intent(in) :: period 
      real(8), intent(in) :: curr_t
      ! local variables
      !integer, parameter :: activation_step = 2
      real(8), parameter :: pi = 3.141592654d0
      
      pos = pos_start
      pos(3) = pos_start(3) + v_mean_up*curr_t + u_amp*sin(curr_t*2*pi/period)

      end subroutine get_source_position_vibration
      !---------------------------------------------------------------
      
      !---------------------------------------------------------------
      subroutine get_Eampl(getEAmpl,noel,npt,ntens,kstep,kinc,time,coords)                                              
         implicit none
         integer, intent(in) :: noel, npt,ntens,kstep,kinc
         real(8),intent(in) :: time(2),coords(3)
         real(8), intent(out) :: getEampl
         real(8), dimension(3) :: dx
         real(8), dimension(3) :: coord_source
         logical :: is_active
         real(8) :: dist ! distance from the external r0 and the current position
         !---- parameters
         real(8), parameter :: r0 = 0.5d0 ! radius within within which epsAmpl is max. and beyond
	                                       ! which starts to decreases
         real(8), parameter :: rs = 0.225d0 !physical radius of the source
         real(8), parameter :: epsAmpl_max = 1.0d-3
         real(8), parameter :: frequency = 30.0d0 !Hz
         real(8), dimension(3) :: head = [0.0d0,0.0d0,3.0d0]   
         !----- end of parameters
		 
		 

         getEAmpl = 0.0d0 ! default value
         call get_source_position(coord_source,is_active,kstep,kinc,time)
         if (.not. is_active) return
		 
		 
         dx = coords - coord_source
         dist = sqrt( dot_product(dx,dx) )
         ! consider a vibratory device with a moving head headheight=3 m height
         ! coord_source is the location of the tip
         dx = coords - coord_source - head/2.0d0
         dx(3) = dx(3)/head(3) 
         dist = sqrt( dot_product(dx,dx) )
         
         ! if(dist< rs) return  ! epsAmpl = 0 within the source
         ! getEampl = min( epsAmpl_max,epsAmpl_max*exp(-1.77d0*(dist-r0)))
         getEampl = epsAmpl_max*min(1.0d0,exp(-1.77d0*(dist-r0)))
		 

      end subroutine get_EAmpl 
      !---------------------------------------------------------------
    