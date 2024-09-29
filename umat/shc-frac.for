! Copyright (C)  2024  Carlos Grandas
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.


! ###  niemunis_tools_nano: A Tensor Operations Library  
! This model relies on a proprietary tensor operations library developed 
! by Andrzej Niemuni at KIT, Karlsruhe. To use this library, please 
! contact Andrzej.Niemunis@kit.edu for access and licensing information.

      !-----------------------------------------------------------------
      module shc_frac
      use umat_tools
      implicit none
      
      
      type MATERIALCONSTANTS
      character(3)::name
      real(8)::
     &    phic,       ! friction angle (Rad)  
     &    A,          ! constant for stiffness      
     &    nu,         ! poisson's ratio
     &    patm,       ! atmospheric pressure
     &    nn,         !  exponent  
     &    Kwater,     ! bulk modulus of water < 2 MPa
     &    CAmpl,      !
     &    CN1,
     &    CN2,
     &    CN3,
     &    Ce,
     &    eref,       ! reference void ratio
     &    Cp,
     &    Cy,
     &    pref,
     &    eampRef,
     &    gSOM,
     &    RSOM
      end type MATERIALCONSTANTS

      type STATEVARIABLES
      real(8)::
     &  epor,
     &  pwater,
     &  gA
          ! other state variables
      
      real(8), dimension(3,3)::     
     &  mb      ! unit tensor perpendicular to the bounding surface B(T)=0
      real(8), dimension(3,3,3,3)::
     &  EE
      end type STATEVARIABLES      
      
       
      private   ! keep it empty: all subroutines/data types will be private to the module by default. 
      public constitutiveRates,get_nQ_nasv,statev2Q,asv2statev,Q2statev,
     &   get_Eampl2,get_dNdt,stiffness,props2mat,MATERIALCONSTANTS,
     &   STATEVARIABLES,get_DepsAccDt,get_mb
      
      !-----------------------------------------------------------------
      contains
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine get_nQ_nasv(nQ,nasv,ntens,nstatv)
      implicit none
      integer,intent(in) :: ntens,nstatv
      integer,intent(out) :: nQ,nasv
      nQ =  2 + 1  ! epor + pwater + gA 
      nasv = 10 ! number of variables for ouput only
      end subroutine get_nQ_nasv
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine statev2Q(statev,nstatv,Q,nQ,ntens)
      implicit none
      integer,intent(in) :: nstatv,nQ,ntens
      real(8),intent(in) :: statev(nstatv)
      real(8),intent(out) :: Q(nQ)
      Q(1:2) = statev(1:2)    ! epor,pore pressure
      Q(3) = statev(11)    ! gA
      end subroutine statev2Q 
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine Q2statev(Q,nQ,statev,nstatv,ntens)
      implicit none
      integer,intent(in) :: nstatv,nQ,ntens
      real(8),intent(out) :: statev(nstatv)
      real(8),intent(in) :: Q(nQ)
      statev(1:2) = Q(1:2)     ! epor,pore pressure
      statev(11) = Q(3)   ! gA
      end subroutine Q2statev 
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine asv2statev(asv,nasv,statev,nstatv)
      implicit none
      integer,intent(in) :: nstatv,nasv
      real(8),intent(out) :: statev(nstatv)
      real(8),intent(in) :: asv(nasv)      
      statev(21:20 + nasv) = asv(:)
      end subroutine asv2statev 
      !-----------------------------------------------------------------
      

      !-----------------------------------------------------------------
      subroutine constitutiveRates(Ttot,Q,D,dt,props,ntens,nQ,nprops,
     &      dTdt,dQdt,asv,nasv,errCode,flg)     
      use niemunis_tools_nano
      implicit none
      integer,intent(in) :: nQ,nprops,nasv,flg,ntens
      real(8),intent(in) :: Ttot(3,3),Q(nQ),D(3,3),props(nprops),dt
      real(8),intent(out) :: dTdt(3,3),dQdt(nQ)
      real(8),intent(inout) :: asv(nasv)
      integer,intent(out) :: errCode
      type(MATERIALCONSTANTS) :: mat
      type(STATEVARIABLES) :: sv	
      real(8) :: T(3,3), DepsAccDt(3,3),EE(3,3,3,3),dgADt,epsAmpl,dNdt  
    
      errCode = 0 ! Everything is ok
      call props2mat('HCA',nprops,props,mat) 
      sv%epor = Q(1)
      sv%pwater = Q(2)
      sv%gA = Q(3)
      epsAmpl = asv(2)
      dNdt = asv(3)          
      T = Ttot + sv%pwater*delta   
      call checkstate(T,sv,mat,errCode)      
      if(errCode .ne. 0) return 
      
      call stiffness(T,mat%A,mat%patm,mat%nn,mat%nu,EE)
      call get_DepsAccDt(T,props,sv%epor,sv%gA,epsAmpl,dNdt,
     &        DepsAccDt,dgADt)
      
      dTdt = (EE .xx. (D - DepsAccDt ) )  


      dQdt(1) = (1.0d0 + sv%epor)*tr(D)  ! d epor/d t
      dQdt(2) = -((1.0d0 + sv%epor)/sv%epor)*mat%Kwater*tr(D)  ! d pwater /d t       
      dQdt(3) = dgAdt  ! hystoriotropy rate
      
      dTdt = dTdt - dQdt(2)*delta   ! total stress rate
      
      return      
      contains
      
      
      

      

       
      !--------------------------------------------------------------------- 
      subroutine checkstate(T,sv,mat,errCode)                                           
      implicit none
      real(8),intent(in) :: T(3,3)
      type(STATEVARIABLES),intent(inout) :: sv
      type(MATERIALCONSTANTS),intent(in) :: mat
      integer,intent(inout) :: errCode
      real(8) :: trT
      
      errCode = 1 ! Something is wrong 
      trT = tr(T)      
      !if(trT>0) call raise('ERROR checkstate: tr(T)>0')      
      if(trT>0) return       
      errCode = 0 ! everything is ok
      end subroutine checkstate
      !--------------------------------------------------------------------- 
      
       
	 




      end subroutine constitutiveRates
      !-----------------------------------------------------------------
      
      
      !--------------------------------------------------------------------- 
      subroutine props2mat(cmname,nprops,props,mat)                     !  Transfers the material properties
      use niemunis_tools_nano                                             !  from vector props to structure mat
      implicit none
      type (MATERIALCONSTANTS), intent (out) :: mat
      integer, intent(in) :: nprops
      real (8), intent(in) :: props(nprops)      
      character(len=*)  :: cmname 
      real(8) :: s,sm 
      
      !mat%name = to_upper(cmname(1:3)) !    make uppercase 
      call to_upper(cmname(1:3),mat%name)  !    make uppercase        
      if(mat%name .ne. 'HCA') call report('ERROR props2mat: '
     &        //'unknown material')  ! stops the calculation  
      !if(nprops<17)  call report('ERROR props2mat: too few props')                                                                    !  ---------HYPOPLASTICITY--------
      if(nprops<15)  call report('ERROR props2mat: too few props')                                                                    !  ---------HYPOPLASTICITY--------
        mat%phic   = props(1)       ! friction angle (Rad)
        mat%A     = props(2)       ! poisson's ratio for increased shear stiffness
        mat%nu     = props(3)       !  granular hardness
        mat%patm     = props(4)       !  exponent in Bauer's compression rule
        mat%nn    = props(5)       !  min. void ratio at p=0
        mat%Kwater    = props(6)       !  crit. void ratio at p=0
        mat%CAmpl    = props(7)       !  max void ratio at p=0 (isotropic compresison)
        mat%CN1  = props(8)       !  baro-pykno exponent for  f_d
        mat%CN2   = props(9)        ! baro-pykno exponent for f_e
        mat%CN3    = props(10)       ! ss-stiffness multiplier
        mat%Ce    = props(11)       ! ss-stiffness multiplier
        mat%eref   = props(12)       ! ss-striffness elastic strain range
        mat%Cp  = props(13)       ! ss-stiffness parameter
        mat%Cy    = props(14)       ! ss-stiffness parameter
        mat%pref = props(15)       ! bulk modulus of water < 2 MPa
        
        ! not used here
        !mat%gSOM = props(16)       ! fraction of initial g_A remaining after a monot. strainig of length R_SOM     
        !mat%RSOM = props(17)  ! not used here
        
      ! DERIVED  
        mat%eampRef = 1.0d-4 ! parameter                  
        
        
      ! CHECK MATERIAL PARAMETERS
      call checkRange(mat%phic,0.0d0,pi/4.0d0,'phic')
      call checkRange(mat%nu,0.0d0,0.5d0,'nu')      
	
      if(mat%phic <= 0 .or. mat%phic*180.0/pi > 85) then
       call report('ERROR props2mat: phic out of range')  ! stops the calculation 
      end if      
      if(any(props<0.0d0)) call report('ERROR props2mat: wrong props' )
      end  subroutine props2mat
      !----------------------------------------------------------------------- 
      
      
      !--------{\large   get\_Fm }---------------------------------------
      subroutine get_Fm(T3in,Fm)
                                                   !  Returns the function $F_M =M/M_c$  by Wolffersdorff,
      use niemunis_tools_nano                                             !  habil-BO p.43. Argument T3in is not normalized
      implicit none
      real (8),intent(in),dimension (1:3,1:3) :: T3in
      real (8),intent(out) :: Fm	  
      real (8) :: t3d(3,3), cos3th, tr2,tr3,tanbe, term1,term2
      T3d=dev(hated(T3in))                                              !   $ \hTb^* $
      tr2 = T3d .xx.T3d                                                 !   $ \hTb^*:\hTb^* $
      tr3 = 3*det(T3d)
      cos3th = 1.0d0                                                    !   If $ \hTb^*:\hTb^* = 0$ No deviatoric stress
      if (tr2.gt.1.d-10)    cos3th = -sq6*tr3/tr2**1.5d0
      if (cos3th.gt.1)  cos3th =  1.0d0
      if (cos3th.lt.-1) cos3th = -1.0d0                                 !   $ \cos 3\theta =-\sqrt{6}\Frac{\tr (\hTb^*\cdot\hTb^*\cdot\hTb^*)}{ [ \hTb^*:\hTb^* ]^{3/2}} $
      tanbe  = sq3*NORM(t3d)
      term1 =  2.0d0+sq2*tanbe*cos3th                                   !  $ \tan\psi =\sqrt{3}\|\hat{\Tb}^*\| $
      term2 =  2.0d0 - tanbe*tanbe
      Fm= dsqrt(abs(tanbe*tanbe/8+term2/term1))-sq2*tanbe*0.25d0    !   $F_M = \sqrt{\frac{1}{8}\tan^2\psi+\Frac{2-\tan^2\psi}{2+\sqrt{2}\tan\psi\cos 3\theta}}-\D\frac{1}{2\sqrt{2}}\tan\psi $
      end subroutine get_Fm   
      
      
      !  --------{\large   stiffness }---------------------------------------      
      subroutine stiffness(Tb,A,patm,nn,nu,LL)
      use niemunis_tools_nano
      implicit none
      real(8), intent(in) :: A,patm,nn,nu
      real(8), intent(inout),dimension(3,3,3,3) :: LL      
      real(8), intent(in) :: Tb(3,3)
      real(8) :: lambda,mu,lambda0,mu0,p,B
      
      p = -tr(Tb)/3.0d0
      lambda0 = A*nu*patm/((1.0d0+nu)*(1.0d0-2.0d0*nu))
      lambda = lambda0*(p/patm)**nn
      mu0 = A*patm/(2.0d0*(1.0d0+nu))
      mu = mu0*(p/patm)**nn      
      LL = lambda*(delta .out. delta) + 2.0d0*mu*Idelta
      
!      B = (-nn/(3.0d0*patm))*(p/patm)**(nn-1.0d0)
!      dLLdT = B*lambda0*((delta .out. delta) .out. delta) + 
!     &        2.0d0*B*mu0*( Idelta .out.  delta)

      end subroutine stiffness       
         




!  --------{\large   get\_mb }---------------------------------------
      subroutine get_mb(Tb,phic,mb)
                                            
      use niemunis_tools_nano                                              
      implicit none
      real (8),intent(in),dimension (1:3,1:3) :: Tb
      real (8),intent(in):: phic      
      real (8),intent(out),dimension (1:3,1:3) :: mb      
      real (8) :: p, q, Mc,  M,Fm,eta,aux
      
      mb = 0
      p = -tr(Tb)/3.0d0
      if (p<=0) return
      q = sqrt(3.0d0/2.0d0)* norm(dev(Tb))
      Mc = 6.0d0*dsin(phic)/(3.0d0-dsin(phic)) 
      call get_Fm(Tb,aux)
      Fm = max(aux ,5.0d0/7.0d0     )      
      M = Mc*Fm                  
      mb = -(p-q*q/(M*M*p))*delta/3.0d0 + 3.0d0*dev(Tb)/(M*M) ! the minus sign in the first term does not coincide with Wichtman et al 2010, CanGeoJour, vol 45,Nr7 eq 4
      !mb = M*M*(pB-2.0d0*p)*delta/3.0d0 +3.0d0*dev(Tb)
      !mb = (1.0d0-q*q/(M*M*p))*delta/3.0d0        
      mb = mb/norm(mb)

      end subroutine get_mb
      
      
!\com ........private function {\large   GET\_FY }.......................................................
      subroutine get_fY(phic,CY,T,fY)
                                                ! \com returns the value of $f_Y$ multiplier that      
      use niemunis_tools_nano
      implicit none      
      real(8),dimension(1:3,1:3),intent(in) :: T
      real(8),intent(in) :: CY,phic
      real(8),intent(out) :: fY
      real(8),dimension(1:3,1:3) :: Td
      real(8) :: I1,I2,I3,Ybar,Y,Yc
      call IIIinvariants(T,I1,I2,I3)                                    !  \com basic invariants of a tensor
      Td =  diagonalized(T)
      Ybar=0                            !  prowizorycznie 12.8.2004
      if(abs(I3)>1.d-15 .and. max(Td(1,1),Td(2,2),Td(3,3)) < 0) then    !  \com usual case: no tension no zero stress
        Y=-I1*I2/I3                                                     !  \com for the Matsuoka -Nakai function
        Yc= (9.0d0-dsin(phic)*dsin(phic))
     &      /(1.0d0-dsin(phic)*dsin(phic))                !  \com for critical state $ Y=Y_c $
        if(Y <= Yc .and. Y >9 ) Ybar=(Y-9.0d0)/(Yc-9.0d0)               !  \com  on the p-axis $ Y=9 $
        if(Y>Yc)  Ybar=1
      endif                                                             ! \com $ \bar{Y} = \Frac{Y-9}{Y_c -9};  Y = -\Frac{I_1 I_2}{I_3};   Y_c = \Frac{9-\sin^2\varphi}{1 - \sin^2\varphi}  $
      !get_fY = mat%CY1*Ybar + Exp(mat%CY2*Ybar)                         ! \com $f_Y  =  \exp(C_4 \bar{Y}^2 ) \mwith C_4 \approx 2 $
      fY =  Exp( CY*Ybar)                         ! \com $f_Y  =  \exp(C_4 \bar{Y}^2 ) \mwith C_4 \approx 2 $
      return                                                            ! \com {\blue physical limit $\epsilon_\acc = \epsilon_\ampl$ will be implemented }
      end subroutine get_fY  
      
      
      
      
      
      ! returns the accumulation strain rate tensor and the updated for a given stress T, void 
      ! ratio, g, and Load frequency N
      subroutine get_DepsAccDt(T,props,epor,gA,epsAmpl,dNdt,
     & DepsAccDt,dgADt)
      use niemunis_tools_nano
      implicit none
      real(8),intent(in) :: T(3,3),epor,dNdt,epsAmpl,props(15),gA
      real(8),intent(out) :: DepsAccDt(3,3),dgAdt
      real(8), parameter ::     pliq= 0.01d0
      real(8) :: phic,A,nu,patm,nn,Kwater,Campl,CN1,CN2,CN3,Ce,gA0,
     &    eref,Cp,Cy,pref,eampRef,s,Yc,p,fAmpl,fNdot,fe,fp,fY,
     &    epsAccdot,fNAdot,fNBdot
      real(8) :: mb(3,3) 
      
      !------- Material constants
      phic   = props(1)
      A      = props(2)
      nu     = props(3)
      patm   = props(4)
      nn     = props(5)
      Kwater = props(6)
      CAmpl  = props(7)
      !CN1    = props(8)*1.2d0    !*1.48d0
      !CN2    = props(9)*1.4d0     ! *1.72d0
      CN1    = props(8)!* 1.48d0
      CN2    = props(9)! *1.72d0
      CN3    = props(10) 
      Ce     = props(11)
      eref   = props(12)
      Cp     = props(13)
      Cy     = props(14)
      pref   = props(15)
      
      eampRef = 1.0d-4 ! parameter
      s = sin(phic)
      Yc = (9.0d0 - s**2)/(1.0d0- s**2)   

      DepsAccDt(:,:) = 0
      dgAdt = 0
      
      p = -tr(T)/3.0d0        
      if(p<0.01d0) return       !stop 'ERROR:HCA p<pliq'
      fAmpl = min( ( epsAmpl/eampRef )**CAmpl, 10.0d0**Campl )
      if(fAmpl >= 100.0d0) fAmpl =100.0d0        
      fNdot = 0.0d0; fNAdot = 0 ; fNBdot = 0
      if (fAmpl > 1.0d-7 ) then ! avoid too small strain amplitudes
          ! the following two lines are for fractal equations
          !CN1 = CN1* 1.48d0
          !CN2 = CN2 *1.72d0
          !gA0 = max(gA,CN1*fAmpl*(2.0d0)**CN2) ! the number in parenthesis was originaly 2.0, why???
          !fNAdot =CN1*CN2*(gA0/(fampl*CN1))** ((CN2-1.0d0)/CN2) 
          !fNdot = fNAdot
          
          ! original formulation
          fNAdot = CN1*CN2*dexp(-gA/(fampl*CN1)) 
          fNBdot = CN1*CN3
          fNdot = fNAdot + fNBdot    
      end if        
      !if(abs(Ce -eref)<1.0d-5) stop 'eref = Ce'
	  if(abs(Ce -eref)<1.0d-5) call report('Error eref = Ce')
      fe = (1.0d0+eref)*(Ce-epor)**2
      fe = fe /( (1.0d0+epor)*(Ce-eref)**2 )        
      fp =dexp(-Cp*(p/pref - 1.0d0)) 
      call get_fY(phic,CY,T,fY)
        
      epsAccdot = fAmpl*fNdot*fp*fe*fY
      !Begin  correction for small stresses     04.08.2014 stellios
      !if (p>0.1d0) then   		  
      !    epsAccdot=epsAccdot*(1.0d0-dexp(-(p-0.1d0)/(2.0d0*0.1d0)))
      !else
      !   epsAccdot=0.0d0;	fNAdot = 0.0d0   
      !endif
      !   Ende  correction for small stresses     04.08.2014      
      call get_mb(T, phic,mb)     
      DepsAccDt = epsAccdot*mb*dNdt                  
      dgAdt = fAmpl*fNAdot *dNdt
      if( dgAdt<0) call report('Error: dgAdt<0')
      
      end subroutine get_DepsAccDt
      
      !------------------------------------------------------------------------
      !SUBROUTINE get_Eampl(Eampl,KSTEP,KINC,TIME,NOEL,NPT,COORDS)      
      !implicit none
      !integer, intent(in) :: kstep,kinc,noel,npt
      !real(8), intent(in) :: time(2),coords(3)
      !real(8), intent(out):: Eampl
      !Eampl = 0
      !Eampl = 1.0d-4   
      !if(kstep == 1) Eampl = Sqrt(3.0d0/2.0d0)*6.0d-4 ! Fig. 3.17 Habil. Torsten
      !if(kstep == 2) Eampl = 1.1d-4
      !if(kstep == 3) Eampl = 1.7d-4
      !if(kstep == 4) Eampl = 2.5d-4
      !end SUBROUTINE get_Eampl
      !------------------------------------------------------------------------
      
      
      !---------------------------------------------------------------------
      subroutine getPile(kstep,fromStep, pileBottom,pileTop,isPile)
      implicit none
      integer,intent(in) :: kstep,fromStep
      real(8),intent(out) :: pileTop(3),pileBottom(3)
      logical,intent(out) :: isPile
      real(8) :: minx,miny,maxx,maxy,c1,c2,cTop,cBottom      
      integer,parameter :: npilesx = 1, npilesy =3
      real(8) :: pilesTop(npilesx,npilesy,3),      
     & pilesBottom(npilesx,npilesy,3)
      integer :: order(fromStep:fromStep +npilesx*npilesy -1,2) ! order(k1,k2) k1th step and the k2th coordinate of the raster
      integer :: i,j,istep,toStep
      
      isPile = .True.
      tostep= fromStep + npilesx*npilesy -1
      if (kstep< fromStep .or. kstep> toStep) then
         isPile = .False.; return
      end if
      
      cTop = -18.6d0; cBottom = -35.0d0  
      minx = 3.0; miny=2.0
      maxx = 3.0; maxy =8.4
      
      istep = fromStep
      do i=1,npilesx
        c1 = minx+ (maxx-minx)*((i-1)/float(npilesx-1))
        do j=1,npilesy
           c2 = miny+ (maxy-miny)*((j-1)/float(npilesy-1))
           pilesTop(i,j,:) =    [c2,cTop,c1]     ! c1,c2,cTop can be put in different positions according to the coordinate system
           pilesBottom(i,j,:) = [c2,cBottom,c1]  ! here it is assumed that z points in the vertical direction of the model.
           order(istep,:) = [i,j]
           istep = istep + 1
        end do
      end do
      
      !write(*,*) 'i, j, pilesCoord'
      !do i=1,npilesx
      !  do j=1,npilesy
      !    write(*,'(i2,i2,3g9.2)') i,j,pilesTop(i,j,:)
      !    write(*,'(i2,i2,3g9.2)') i,j,pilesBottom(i,j,:)
      !  end do
      !end do
      
       ! the array Order can be changed to install piles in   different sequences 
      order(fromStep,:)    = [1,1]
      order(fromStep+1,:)  = [1,2]      
      order(fromStep+2,:)  = [1,3]
            	  
      i = order(kstep ,1) ;j = order(kstep ,2)  
      pileBottom  = pilesBottom(i,j,:)
      pileTop  = pilesTop(i,j,:)
      end subroutine getPile 
      !---------------------------------------------------------------------

	 


      !--------{\large   getEampl}---------------------------------------    
      subroutine get_Eampl2(getEampl,noel,npt,ntens,kstep,kinc,time,
     &    coords)                                              
      use niemunis_tools_nano                                      
      implicit none
      integer, intent(in) :: noel, npt,ntens,kstep,kinc
      real(8),intent(in) :: time(2),coords(3)
      real(8),intent(out) :: getEampl
      real(8) :: xy_pile_top(3),xy_pile_bottom(3),xy_pile_current(3),
     &   r0,epsAmpl_max,frequency,penetrationRate,penetrationLength,
     &   dL,alpha,dx,lfac,rp
      integer:: pileInstallStep      
      logical :: isPile

      pileInstallStep = 8 ! in this step the pile will be vibrated into the soil
      getEampl = 0 ! default value for all times where no accumulation occurs
      !if( kstep>1) write(*,*) 'before getPile. kstep,kinc',kstep,kinc
      call getPile(kstep,pileInstallStep,xy_pile_bottom,xy_pile_top,
     &        isPile)
      !      if(lfac < 0) return  ! return getEAmpl = 0 for steps different than 11,12,13
       !if( kstep>1) write(*,*) 'after getPile. isPile=', isPile
       if (.not. isPile) return
	  
!      if(lfac < 0) return  ! return getEAmpl = 0 for steps different than 11,12,13
!      if(kstep .ne. pileInstallStep) return 	  
	  
	  ! The strain amplitude near to the pile tip is 1.0d-3 and decreases
	  ! exponentially with the distance. 
	  
	  !xy_pile_top(:) = [2.0d0,-18.6d0]
	  !xy_pile_bottom(:) = [2.0d0,-31.6d0]
      r0 = 0.5d0 ! radius within within which epsAmpl is max. and beyond
	             ! which starts to decreases
      rp = 0.15d0 !pile radius
      epsAmpl_max = 1.0d-3
      frequency = 34.0d0 !Hz
      penetrationRate = 1.0d0 ! m/s
      penetrationLength= sqrt( dot_product((xy_pile_top-xy_pile_bottom),
     &                       (xy_pile_top-xy_pile_bottom)) )
	  dL = time(1)*penetrationRate 
      !if( kstep>1) write(*,*) 'getEampl before dL<penetra...'      
      if(dL < penetrationLength) then
        xy_pile_current = xy_pile_top ! moving tip
     &	              + dL*(xy_pile_bottom-xy_pile_top)/penetrationLength	         		 
        	
        !xy_pile_current = xy_pile_top+(xy_pile_bottom-xy_pile_top)*lfac        
	    ! distance between the current gauss point and the current position
		! of the pile tip
        dx = sqrt( dot_product((xy_pile_current-coords),  
     &                         (xy_pile_current-coords)) )	 
        if(dx< rp) return  ! epsAmpl = 0 within the pile
        getEampl = min( epsAmpl_max,epsAmpl_max*exp(-1.77d0*(dx-r0)))
      end if
      !if( kstep>1) write(*,*) 'end getEampl=', getEampl
	  
	  !*********************
      ! getEampl = epsAmpl_max
	  !*********************
	  
      end subroutine get_Eampl2
      !------------------------------------------------------------------------
   
      
      !------------------------------------------------------------------------
      SUBROUTINE get_dNdt(dNdt,KSTEP,KINC,TIME,NOEL,NPT,COORDS)
      implicit none
      integer, intent(in) :: kstep,kinc,noel,npt
      real(8), intent(in) :: time(2),coords(3)
      real(8), intent(out):: dNdt
      dNdt = 30.0d0     
      end SUBROUTINE get_dNdt
      !------------------------------------------------------------------------
      
       
      
      
      end module shc_frac
      !-----------------------------------------------------------------