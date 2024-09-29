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
      module wo2_sch
      use umat_tools
      implicit none
      
      type MATERIALCONSTANTS
      character(3)::name
      real(8)::
     &    phic,       ! friction angle (Rad)        
     &    nu,         ! poisson's ratio
     &    hs,         !  granular hardness
     &    en,         !  exponent in Bauer's compression rule
     &    ed0,        !  min. void ratio at p=0
     &    ec0,        !  crit. void ratio at p=0
     &    ei0,        !  max void ratio at p=0 (isotropic compresison)
     &    alpha,      !  baro-pykno exponent for  f_d
     &    beta,       ! baro-pykno exponent for f_e
     &    m_2,        ! ss-stiffness multiplier
     &    m_5,        ! ss-stiffness multiplier
     &    Rmax,       ! ss-striffness elastic strain range
     &    betax,      ! ss-stiffness parameter
     &    Chi,        ! ss-stiffness parameter
     &    Kwater,     ! bulk modulus of water < 2 MPa
     &    epor0,      ! void ratio  if trT = 0 
     &    etaBag,     ! Bagnold's viscosity constant
     &    KBag,       ! Bagnold's viscosity constant
     
     &    prefsmallE, ! reference pressure for small stiffness
     &    xismallE,   ! exponent for small stiffness
     &    nsmallE,    ! barotropy exponent for small stiffness
     
     
          ! other mat constants
     &    pt,     
     &    phimax,
     &    Mcmax,  ! maximum value of F_B (a function of phimax)
     &    az,
     &    az2,
     &    Mc,
     &    Me,
     &    fi,
     &    b2,
     &    hi
      end type MATERIALCONSTANTS

      type STATEVARIABLES
      real(8)::
     &  epor,
     &  pwater,
     &  Bauer,
     &  ed,
     &  ec,
     &  ei,
          ! other state variables
     &  Fm,      
     &  fb,
     &  fe,
     &  fs,
     &  fd,    
     &  p,
     &  q,
     &  B,
     &  M,
     &  rho,  ! mobilization of intergranular strain    
     &  normBB,
     &  Y 
!     &  gA
      real(8), dimension(3,3)::
     &  NN,       
     &  h,       ! intergranular strain
     &  hunit,
!     &  dBdT,    ! tensor perpendicular to the bounding surface B(T)=0
!     &  nb,      ! unit tensor perpendicular to the bounding surface B(T)=0
     &  Tvis,     ! Viscous stress: Ttot = Teff - poreWater + Tvis
     &  NNhat,
     &  BB,
     &  mb
      real(8), dimension(3,3,3,3)::
     &  LLhat,      
     &  LL,     
     &  MM,
     &  LLhat_inv
      end type STATEVARIABLES      
      
       
      private   ! keep it empty: all subroutines/data types will be private to the module by default. 
      public constitutiveRates,get_nQ_nasv,statev2Q,asv2statev,Q2statev,
     &   statev2asv 
      
      !-----------------------------------------------------------------
      contains
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine get_nQ_nasv(nQ,nasv,ntens,nstatv)
      implicit none
      integer,intent(in) :: ntens,nstatv
      integer,intent(out) :: nQ,nasv
      nQ =  2 + ntens +1 +1 +ntens ! epor + pwater + ntens of intergranular strain + gA +epsAcc
      nasv = 20 ! number of variables for ouput only
      end subroutine get_nQ_nasv
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine statev2Q(statev,nstatv,Q,nQ,ntens)
      implicit none
      integer,intent(in) :: nstatv,nQ,ntens
      real(8),intent(in) :: statev(nstatv)
      real(8),intent(out) :: Q(nQ)
      Q(1:2) = statev(1:2)    ! epor,pore pressure
      Q(3:2+ntens) = statev(3:2+ntens)    ! intergranular strain
      Q(3+ntens) = statev(11) !gA
      Q(4+ntens) = statev(12) !epsAcc
      Q(5+ntens:4+2*ntens) = statev(13:12+ntens)    ! viscous stress
      end subroutine statev2Q 
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine Q2statev(Q,nQ,statev,nstatv,ntens)
      implicit none
      integer,intent(in) :: nstatv,nQ,ntens
      real(8),intent(out) :: statev(nstatv)
      real(8),intent(in) :: Q(nQ)
      statev(1:2) = Q(1:2)     ! epor,pore pressure
      statev(3:2+ntens) = Q(3:2+ntens)   ! intergranular strain
      statev(11) = Q(3+ntens) !gA
      statev(12) = Q(4+ntens) !epsAcc
      statev(13:12+ntens) = Q(5+ntens:4+2*ntens)    ! viscous stress
      end subroutine Q2statev 
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine  statev2asv(statev,nstatv,asv,nasv)
      implicit none
      integer,intent(in) :: nstatv,nasv
      real(8),intent(in) :: statev(nstatv)
      real(8),intent(out) :: asv(nasv)      
      asv(:) =statev(21:20 + nasv)
      end subroutine statev2asv 
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
      use umat_tools ! rkinc, rkstep,....
      use shc_frac, stiffness_HCA => stiffness,props2mat_HCA=>props2mat,
     &    MATERIALCONSTANTS_HCA => MATERIALCONSTANTS,
     &    STATEVARIABLES_HCA => STATEVARIABLES,
     &    get_mb_HCA => get_mb
      implicit none
      integer,intent(in) :: nQ,nprops,nasv,flg,ntens
      real(8),intent(in) :: Ttot(3,3),Q(nQ),D(3,3),props(nprops),dt
      real(8),intent(out) :: dTdt(3,3),dQdt(nQ),asv(nasv)
      integer,intent(out) :: errCode
      type(MATERIALCONSTANTS) :: mat
      type(STATEVARIABLES) :: sv
      type(MATERIALCONSTANTS_HCA) :: mat_h
      type(STATEVARIABLES_HCA) :: sv_h    
      real(8) :: T(3,3), hD,A,B,C,rhox,dhdt(3,3),huhu(3,3,3,3),hnorm,
     &    EE_HCA(3,3,3,3),AAA(3,3),BBB(3,3),CCC(3,3),Dh(3,3),xi
      real(8) :: epsAmpl,dNdt, DepsAccDt(3,3),dgADt,EE(3,3,3,3), 
     &  dgAdt_mon,aux,mm(3,3),checkE(3,3,3,3),Tvisn(3,3),dTvisdt(3,3),
     &  dpwdt,depordt
      integer :: nprops_h
      real(8),allocatable :: props_h
      integer,parameter :: solution =1, jacobian =2,jacobian0 =3
            
      if (flg==jacobian) then
        continue
      end if
    
      errCode = 0 ! Everything is ok
      call props2mat('WO2',nprops,props,mat)       
      sv%epor = Q(1)
      sv%pwater = Q(2)
      sv%h = map2D(Q(3:2+ntens),ntens)  ! intergranular strain (4-6 components)
      sv_h%gA = Q(3+ntens)  ! gA  
      !epsAmpl = asv(8)
      !dNdt = asv(9)
      Tvisn = map2T(Q(5+ntens:4+2*ntens),ntens) ! viscous stress (4-6 components)
      
    
      T = Ttot + sv%pwater*delta !- Tvisn  
      call checkstate(T,sv,mat,errCode)      
      if(errCode .ne. 0) return      
      
      
      !****** from HCA********
      call props2mat_hca('HCA',nprops-16,props(17:),mat_h)
      call stiffness_hca(T,mat_h%A,mat_h%patm,mat_h%nn,mat_h%nu,EE_hca)      
      !call get_Eampl(epsAmpl,rKSTEP,rKINC,rTIME,rNOEL,rNPT,rCOORDS) ! all this variables are passed in umat_tools module
      
      call get_Eampl(epsAmpl,rNOEL,rNPT,ntens,rKSTEP,rKINC,rTIME,
     &   rCOORDS) ! get pile... many piles in 3 D Postdamer Platz
      
      call get_dNdt(dNdt,rKSTEP,rKINC,rTIME,rNOEL,rNPT,rCOORDS)      
      call get_DepsAccDt(T,props(17:),sv%epor,sv_h%gA,epsAmpl,dNdt,
     &        DepsAccDt,dgADt)
      call get_mb_HCA(T, mat_h%phic,sv_h%mb)
      !***** end from HCA
     
      Dh = D - DepsAccDt
      xi = 1.0d0 ! monotonic behaviour      
      if (norm(Dh) +norm(D)>0) xi = (1.0d0- 
     &       norm(DepsAccDt)/( norm(D)+norm(Dh)) )      
      
	  ! 19.12.2023
	  xi = 1.0d0 ! assume monotonic stiffness and flow rule!
	  ! 19.12.2023
	  
      hNorm = norm(sv%h)
      if(hNorm > mat%Rmax) sv%h = sv%h*mat%Rmax/hNorm  ! control to avoid sv%rho>1 
      sv%rho = norm(sv%h)/mat%Rmax
      if(sv%rho < 1.0d0 - 1.0d-4) errCode = 2  ! use reduced stress tolerance

      call getLLNN(T,mat,sv)      
      if (mat%m_5 < 2 ) then ! ignore intergranular strain effect      
        !dTdt = (sv%LL .xx. D) + sv%NN*norm(D)  ! effective stress rate
        
        call getYmb(T,mat,sv)        
        aux = sv%mb .xx. normalized(DepsAccDt)        
        EE = xi*sv%LL + (1.0d0 -xi)*EE_hca/(1.0d0 + sv%Y)        
        call smallPressureStiffness(T,mat,EE)        
        mm = xi*sv%mb + (1.0d0 -xi)*sv_h%mb 
        dTdt = ( EE .xx. (Dh - sv%Y* mm*norm(Dh) ))        
        
        !dTdt = ( EE_hca .xx. (D   - DepsAccDt ))        
        ! for debug
        !AAA = -sv%Y*(sv%LL .xx. sv%mb)*norm(D) - sv%NN*norm(D)
        !BBB = ( sv%LL .xx. (D - sv%Y*sv%mb*norm(D) ))
        !CCC = dTdt - BBB
        !if(any((CCC>1.0d-10) .or. (CCC<-1.0d-10) ) ) stop 'bad CCC'
        !continue       
      end if
      
      if (mat%m_5 >= 2) then  ! consider intergranular strain effect
        !call report('Error: Interg. strain not implemented yet')
        sv%hunit = normalized(sv%h)        
        hD = (sv%hunit .xx. Dh)
        if(hD > 0) then
           dhdt = ((jdelta - 
     &            (sv%hunit .out. sv%hunit)*sv%rho**mat%betax) .xx. Dh) ! evolution of h with Dh not D
        end if
        if(hD <= 0)    dhdt = Dh  ! reversal      
        sv%rho = norm(sv%h)/mat%Rmax        
        rhox = sv%rho**mat%chi
        A = rhox*mat%m_2 + (1.0d0-rhox)*mat%m_5
        B = rhox*(1.0d0-mat%m_2)
        C = rhox*(mat%m_5-mat%m_2)
        huhu = (sv%hunit .out. sv%hunit)
        sv%MM = A*sv%LL
        if(hD>0) then        
          sv%MM = sv%MM + B*(sv%LL .xx. huhu) 
     &           + rhox* (sv%NN .out. sv%hunit)          
        end if
        if(hD <= 0) then  ! reversal
          sv%MM = sv%MM + C*(sv%LL .xx. huhu)           
        end if
        EE = xi*sv%MM + (1.0d0 -xi)*EE_hca
        call smallPressureStiffness(T,mat,EE) 
        dTdt = (EE .xx. Dh)        
        !dTdt = (sv%MM .xx. D)
      end if ! end consider intergranular strain effect
      
      if( (tr(T) >= 0) .and. tr(D) < 0) dTdt = D ! allow effective stress evolution at p=0 for compressive strain rates 
          
      
      
      depordt = (1.0d0 + sv%epor)*tr(D)                  ! void ratio rate 
      dpwdt = -((1.0d0 + sv%epor)/sv%epor)*mat%Kwater *tr(D) ! pore water pressure rate  
      !dgAdt_mon = (log(mat_h%gSOM)/mat_h%RSOM)*sv_h%gA*norm(Dh)*dNdt   
      !dgAdt_mon = sv_h%gA*6.0d2* norm(Dh- ( Dh .xx. sv_h%mb)*sv_h%mb) ! only Dh perpendicular to mb reduces gA
      !dgAdt_mon = sv_h%gA*2.0d2* norm(Dh)
      !!!! not used here
      !!!!!dgAdt_mon = sv_h%gA*2.0d2* norm(Dh) 
      sv%Tvis = mat%etaBag*(dev(D)-mat%KBag*delta*norm(dev(D))/sq3)
      dTvisdt = 0
      if(dt> 1.0d-8) dTvisdt = (sv%Tvis - Tvisn)/dt    ! viscous stress rate 
      
      dTdt = dTdt - dpwdt*delta !+ dTvisdt ! total stress rate
      
      dQdt(:) = 0            ! rate of change of state variables 
      dQdt(1) = depordt  ! d epor/d t
      dQdt(2) = dpwdt  ! d pwater /d t       
      dQdt(3:2+ntens) = map2stran(dhdt,ntens)  ! intergranular strain  d h / d t                
      dQdt(3+ntens) =    dgAdt  !-   dgAdt_mon
      dQdt(4+ntens) = norm(DepsAccDt) 
      dQdt(5+ntens:4+2*ntens) = map2stress(sv%Tvis,ntens)  ! viscous stress ratio
      
      ! additional info variables. asv(1) is reserved for the number of RK subinc!
      asv(:) = 0      
      asv(2) = sv%rho
      call getPQ(T,sv%p,sv%q)
      asv(4) = sv%p  !
      asv(5) = sv%q  !
      sv%M  = get_Mtheta(T,mat%phic)
      asv(6) = sv%M
      if(abs(sv%p*sv%M)>1.0d-8) asv(7) = sv%q/(sv%p*sv%M ) !
      asv(8) = xi !
      asv(9) = norm(D) 
      asv(10) = norm(Dh) 
      asv(11) = norm(DepsAccDt)
      asv(12) = epsAmpl	  
      
      return      
      contains     
       
      

      !--------------------------------------------------------------------- 
      subroutine props2mat(cmname,nprops,props,mat)                     !  Transfers the material properties
      use niemunis_tools_nano                                             !  from vector props to structure mat
      use umat_tools
      implicit none
      type (MATERIALCONSTANTS), intent (out) :: mat
      integer, intent(in) :: nprops
      real (8), intent(in) :: props(nprops)      
      character(len=*)  :: cmname 
      real(8) :: s,sm
      
      !mat%name = to_upper(cmname(1:3)) !    make uppercase  
      call to_upper(cmname(1:3),mat%name)  !    make uppercase        

      if(mat%name .ne. 'WO2') call report('ERROR props2mat: '
     &        //'unknown material '//trim(cmname))  ! stops the calculation
        !  ---------HYPOPLASTICITY--------
        mat%phic   = props(1)       ! friction angle (Rad)
        mat%nu     = props(2)       ! poisson's ratio for increased shear stiffness
        mat%hs     = props(3)       !  granular hardness
        mat%en     = props(4)       !  exponent in Bauer's compression rule
        mat%ed0    = props(5)       !  min. void ratio at p=0
        mat%ec0    = props(6)       !  crit. void ratio at p=0
        mat%ei0    = props(7)       !  max void ratio at p=0 (isotropic compresison)
        mat%alpha  = props(8)       !  baro-pykno exponent for  f_d
        mat%beta   = props(9)        ! baro-pykno exponent for f_e
        mat%m_2    = props(10)       ! ss-stiffness multiplier
        mat%m_5    = props(11)       ! ss-stiffness multiplier
        mat%Rmax   = props(12)       ! ss-striffness elastic strain range
        mat%betax  = props(13)       ! ss-stiffness parameter
        mat%Chi    = props(14)       ! ss-stiffness parameter
        mat%Kwater = props(15)       ! bulk modulus of water < 2 MPa
        !mat%epor0  = props(16)       ! void ratio  if trT = 0
        
        
        ! small stiffness (constant for all materials)
        mat%xismallE = 2.0d0
        mat%nsmallE = 0.4d0
        mat%prefsmallE = 1.0d0
        
        ! Bagnold viscosity
        mat%etaBag =  0.02d0 ! kPa s 
        mat%KBag = 0.7d0   
        
      ! DERIVED  
                  
        
        s = dsin(mat%phic)
        mat%az = sq3*(3.0d0-s)/(2.0d0*sq2*s)                          !  $\sqrt{3}(3-\sin{\phi_c})/(2 \sqrt{2} \sin{\phi_c})$
        mat%az2 = mat%az*mat%az
        mat%Me = 6.0d0*s/(3.0d0+s)                                    !   $M_E = \dfrac{6 \sin{\phi_c}}{3+\sin{\phi_c}}$
        mat%fi = (9.0d0-s*s)/(1.0d0-s*s)                              !   $\phi = \dfrac{9-\sin^2 \varphi_c}{1-\sin^2 \varphi_c}$
        mat%b2 = (1 + mat%az2/3 + mat%az/sq3)*(1-2*mat%nu)/(1+mat%nu)-1   ! \com  $b^2$ for increased shear stiffness               ! \com  $b^2$ for increased shear stiff
        mat%hi= 3+mat%az*mat%az-mat%az*sq3 *(((mat%ei0-mat%ed0)/          !  \com habil-BO (2.72) $h_i=\left[ 3+a^2-a\sqrt{3}\left(\D\frac{e_{i0}-e_{d0}}{e_{c0}-e_{d0}}\;\right)^{\alpha}\right]^{-1}$
     &          (mat%ec0-mat%ed0))**mat%alpha )
      ! CHECK MATERIAL PARAMETERS
      
      call checkRange(mat%phic,0.0d0,pi/4.0d0,'phic')
      call checkRange(mat%nu,0.0d0,0.5d0,'nu')
      call checkRange(mat%hs,0.0d0,1.0d20,'hs')
      call checkRange(mat%en,0.0d0,1.0d20,'en')
      call checkRange(mat%ed0,0.0d0,1.0d20,'ed0')
      call checkRange(mat%ec0,0.0d0,1.0d20,'ec0')
      call checkRange(mat%ei0,0.0d0,1.0d20,'ei0')
      call checkRange(mat%phic,0.0d0,85.0d0*pi/180.0d0,'phic')

      if(any(props<0.0d0)) call report('ERROR props2mat: wrong props')
      
      end  subroutine props2mat
      !----------------------------------------------------------------------- 

       

      !--------------------------------------------------------------------- 
      subroutine checkstate(T,sv,mat,errCode)    
      !use niemunis_tools_nano                                            
      implicit none
      real(8),intent(in) :: T(3,3)
      type(STATEVARIABLES),intent(inout) :: sv
      type(MATERIALCONSTANTS),intent(in) :: mat
      integer,intent(inout) :: errCode
      real(8) :: trT
      
      errCode = 1 ! Something is wrong 
      trT = tr(T)
      if(trT> -5.0d0 ) then
      continue
      end if
      !if(trT>0) call raise('ERROR checkstate: tr(T)>0')      
      if(trT>0) return    
      sv%Bauer=dexp(-(-trT/mat%hs)**mat%en)  ! this will be calculated in getLLNN again
      sv%ed = mat%ed0*sv%Bauer
      sv%ec = mat%ec0*sv%Bauer
      sv%ei = mat%ei0*sv%Bauer
      !if(sv%epor < sv%ed)   call report('ERROR checkstate: e<ed')
      !if(sv%epor > sv%ei)   
      !call raise('ERROR checkstate: e>ei')      
      errCode = 0 ! everything is ok
      end subroutine checkstate
      !--------------------------------------------------------------------- 
      
      !-----------------------------------------------------------------------
      subroutine getLLNN(Tb,mat,sv)
      use niemunis_tools_nano
      implicit none
      type (MATERIALCONSTANTS),intent(in) :: mat
      type (STATEVARIABLES),intent(inout) :: sv
      real(8), intent(in) :: Tb(3,3)
      real(8) :: trT,That(3,3),aux,hi,fbfefs,LLold(3,3,3,3),mtd,med,
     &          fdq 
      logical :: ok
      
      trT = tr(Tb)
      sv%LL = 0.0d0 ; sv%NN = 0.0d0
      if(trT >= 0) return
      
      sv%Bauer=dexp(-(-trT/mat%hs)**mat%en)  
      sv%ed = mat%ed0*sv%Bauer
      sv%ec = mat%ec0*sv%Bauer
      sv%ei = mat%ei0*sv%Bauer
      
      sv%Fm = get_Fm(Tb)                                              !   $F_M(\Tb)$ limit stress obliquity (depends on $\theta$)
      That = hated(Tb)                                                !   $ \hat {\Tb} = \Tb / \tr \Tb
      sv%LLhat= sv%Fm*sv%Fm*Jdelta+mat%az2*(That .out. That)      !   linear hp stiffness $ \hat{\cE} = a^2 \left[ \left(\Frac{F_M}{a}\right)^2 \cI + \hTb \hTb \right] $    
      
      aux = ((mat%ei0-mat%ed0)/(mat%ec0-mat%ed0))**mat%alpha
      hi= 3.0d0+mat%az2-mat%az*sq3*aux  ! 2.72
      if(hi .le. 0.0d0) call report('ERROR getLLNN: hi<=0')
      !if(hi .le. 0.0d0) stop 'ERROR:getLLNN stop: hi<=0'
      aux = (mat%hs/mat%en)*(mat%ei0/mat%ec0)**mat%beta
      sv%fb=aux*(1.d0+sv%ei)/sv%ei*(-trT/mat%hs)**(1.d0-mat%en)/hi
      sv%fe=(sv%ec/sv%epor)**mat%beta        
      sv%fs= 1.0d0 /(That .xx. That) ! if trT ==0 -> fs=3 avoid null division.
      fbfefs = sv%fb*sv%fe*sv%fs
      LLold= fbfefs*sv%LLhat
      sv%LL = LLold
      
      if (sv%epor .gt. sv%ed)then
        sv%fd = ( (sv%epor-sv%ed)/(sv%ec-sv%ed) )**mat%alpha
      end if
      if (sv%epor .lt. sv%ed) then
        sv%fd = -( (sv%ed-sv%epor)/(sv%ec-sv%ed) )**mat%alpha
      end if            
      mtd = -sv%ed/mat%hs*mat%en*(-trT/mat%hs)**(mat%en-1)              ! \com cf. Habil-BO (4.222) fd   consistent with the lower bound
      med =  1.0d0                                                      !\com $ M_T^{(d)} = \pp{F_d}{\tr\Tb} = -\frac{e_d}{h_s} n \left(\frac{-\tr\Tb}{h_s}\right)^{n-1}\, ;  M_e^{(d)}  = 1$.
      fdq=-(med*sq3*(1+sv%epor)+ mtd*sv%fb*sv%fe*sq3*(3+mat%az2))             !\com $\bar{f}_d =  -\Frac{M_e^{(d)}  \sqrt3~(1 + e)   + M_T^{(d)} f_b f_e \Frac{3}{\sqrt3}(3+a^2)}{ M_T^{(d)} f_b f_e 3a }$
     &      /(mtd*3*sv%fb*sv%fe*mat%az)
      if(sv%fd.lt.1)  sv%fd= sv%fd + (1-sv%fd)**5 * fdq   !\com $f_d = {\left(\frac{e - e_d}{e_c - e_d}\right)}^\alpha  + \left[1 -  {\left(\frac{e - e_d}{e_c - e_d}\right)}^\alpha\right]^z \bar{f}_d \, .$

      sv%NN = fbfefs*sv%fd*mat%az*sv%Fm*( That + dev(That))

      
      ! increased shear stiffness with poisson ratio
      !sv%LL=LLold + fbfefs*mat%b2*(Jdelta - (delta .out. delta)/3.0d0)
      !sv%NN = sv%LL .xx. (Inv(LLold,ok) .xx. sv%NN)  ! modify N to preserve the flow rule direction $\Bb = \cL^{-1} : \Nb$
      
      end subroutine getLLNN 
      !--------------------------------------------------------------------- 
      
      !-----------------------------------------------------------------------
      subroutine getYmb(Tb,mat,sv)
      use niemunis_tools_nano
      implicit none
      type (MATERIALCONSTANTS),intent(in) :: mat
      type (STATEVARIABLES),intent(inout) :: sv
      real(8), intent(in) :: Tb(3,3)
      real(8) :: trT,That(3,3),aux,hi,fbfefs,LLold(3,3,3,3),mtd,med,
     &          fdq, fac
      logical :: ok
      
      trT = tr(Tb)
      sv%Bauer=dexp(-(-trT/mat%hs)**mat%en)  
      sv%ed = mat%ed0*sv%Bauer
      sv%ec = mat%ec0*sv%Bauer
      sv%ei = mat%ei0*sv%Bauer
      
      sv%Fm = get_Fm(Tb)                                              !   $F_M(\Tb)$ limit stress obliquity (depends on $\theta$)
      That = hated(Tb)                                                !   $ \hat {\Tb} = \Tb / \tr \Tb
      sv%LLhat= sv%Fm*sv%Fm*Jdelta+mat%az2*(That .out. That)      !   linear hp stiffness $ \hat{\cE} = a^2 \left[ \left(\Frac{F_M}{a}\right)^2 \cI + \hTb \hTb \right] $    
      fac = (sv%Fm/mat%az)**2 + (That .xx. That)      
      sv%LLhat_inv = (Jdelta-(That .out. That)/fac)/(sv%Fm)**2      
      sv%NNhat = mat%az*sv%Fm*( That + dev(That))
      sv%BB = (sv%LLhat_inv .xx. sv%NNhat)
      sv%normBB = norm(sv%BB)
      sv%mb = -normalized(sv%BB)       
           
      if (sv%epor .gt. sv%ed)then
        sv%fd = ( (sv%epor-sv%ed)/(sv%ec-sv%ed) )**mat%alpha
      end if
      if (sv%epor .lt. sv%ed) then
        sv%fd = -( (sv%ed-sv%epor)/(sv%ec-sv%ed) )**mat%alpha
      end if      
      
      mtd = -sv%ed/mat%hs*mat%en*(-trT/mat%hs)**(mat%en-1)              ! \com cf. Habil-BO (4.222) fd   consistent with the lower bound
      med =  1.0d0                                                      !\com $ M_T^{(d)} = \pp{F_d}{\tr\Tb} = -\frac{e_d}{h_s} n \left(\frac{-\tr\Tb}{h_s}\right)^{n-1}\, ;  M_e^{(d)}  = 1$.
      fdq=-(med*sq3*(1+sv%epor)+ mtd*sv%fb*sv%fe*sq3*(3+mat%az2))             !\com $\bar{f}_d =  -\Frac{M_e^{(d)}  \sqrt3~(1 + e)   + M_T^{(d)} f_b f_e \Frac{3}{\sqrt3}(3+a^2)}{ M_T^{(d)} f_b f_e 3a }$
     &      /(mtd*3*sv%fb*sv%fe*mat%az)
      if(sv%fd.lt.1)  sv%fd= sv%fd + (1-sv%fd)**5 * fdq   !\com $f_d = {\left(\frac{e - e_d}{e_c - e_d}\right)}^\alpha  + \left[1 -  {\left(\frac{e - e_d}{e_c - e_d}\right)}^\alpha\right]^z \bar{f}_d \, .$

      sv%Y = sv%fd*sv%normBB

      
      ! increased shear stiffness with poisson ratio
      !sv%LL=LLold + fbfefs*mat%b2*(Jdelta - (delta .out. delta)/3.0d0)
      !sv%NN = sv%LL .xx. (Inv(LLold,ok) .xx. sv%NN)  ! modify N to preserve the flow rule direction $\Bb = \cL^{-1} : \Nb$
      
      end subroutine getYmb 
      !--------------------------------------------------------------------- 
      
      
      

      !---------------------------------------------------------------------
      function get_Fm(T3in)                                             !  Returns the function $F_M =M/M_c$  by Wolffersdorff,
      !use niemunis_tools_nano                                             !  habil-BO p.43. Argument T3in is not normalized
      implicit none
      real (8),intent(in),dimension (1:3,1:3) :: T3in
      real (8) :: GET_Fm
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
      get_Fm= dsqrt(abs(tanbe*tanbe/8+term2/term1))-sq2*tanbe*0.25d0    !   $F_M = \sqrt{\frac{1}{8}\tan^2\psi+\Frac{2-\tan^2\psi}{2+\sqrt{2}\tan\psi\cos 3\theta}}-\D\frac{1}{2\sqrt{2}}\tan\psi $
      end function get_Fm
      !-----------------------------------------------------------------
      
      
      
      
      !-----------------------------------------------------------------  
      subroutine getPQ(Tb,p,q)
      use niemunis_tools_nano
      implicit none            
      real(8), intent(in) :: Tb(3,3)
      real(8),intent(out) :: p,q            
            
      p = -tr(Tb)/3.0d0
      q = dsqrt(3.0d0/2.0d0)*norm(dev(Tb))
      end subroutine getPQ
      !-----------------------------------------------------------------  
      
       

      !-----------------------------------------------------------------
      function get_Mtheta(T,phi)   ! returns the slope M=q/p of the CSL for a given T and phi      
      use niemunis_tools_nano        ! according to Matsuoka-Nakai criterium: $M(\theta)$ from quadratic or qubic equation
      implicit none
      real(8) :: get_Mtheta                                             !  $M(\theta)$
      real(8), intent(in)::  T(3,3), phi
      real(8) :: c3t,c3t2,cos3theta,fi, fi2, fi3, Hc3t,alpha3,num 

      c3t= -3.0d0*sq6*det(  normalized( dev(T) )    )                   ! $\cos(3 \theta)= - 3 \sqrt6 \det \vec{\Tb}^*$
      c3t =  max(min(c3t,1.0d0),-1.0d0)
      c3t2 = c3t*c3t

      fi = (9.0d0-dsin(phi)**2) /(1.0d0-dsin(phi)**2)
      if ( abs(c3t)< 1.0d-7 ) then                                       ! special square solution:
        cos3theta = c3t        
        get_Mtheta = sqrt(3*(fi-9.0d0)/(fi-3.0d0))              
        return
      end if
                                                                       ! qubic solution:
      fi2 = fi*fi
      fi3 = fi*fi*fi
      if(c3t > 0)  Hc3t = 1.0d0                                         ! Heaviside function $H(\cos(3 \theta))$
      if (c3t < 0) Hc3t = 0.0d0
      num = -27 +27*fi - 9*fi2 + 18*c3t2*fi2 + fi3 - 2*c3t2*fi3         ! numerator of $   \dfrac{ -27 + 27  \phi - 9  \phi^2 + 18 \cos^2(3\theta)  \phi^2 + \phi^3 -  2 \cos^2(3\theta)  \phi^3 } { (-3 + \phi)^3 }    $
      alpha3 = acos(  sign(1.0d0,c3t)*num/(-3.0d0+fi)**3  )/3.0d0       ! $\alpha = \arccos\left\{ \sign\left[\cos(3\theta)\right] \right. $
      get_Mtheta =((3*fi-9)/(fi*c3t))*(0.5d0-cos(alpha3+ pi*Hc3t/3))    ! $M =  \frac{3 \phi  - 9  }{ \phi  \cos(3\theta) }  \left[\frac12 -\cos\left[ \frac{\alpha}3 + \frac{\pi}3  H( \cos(3\theta) ) \right]  \right]    $      
      end function get_Mtheta
      !-----------------------------------------------------------------
      
    
      !-----------------------------------------------------------------  
      subroutine smallpressurestiffness(T,mat,EE) ! after osinov
      use niemunis_tools_nano
      implicit none            
      type (MATERIALCONSTANTS),intent(in) :: mat
      real(8), intent(in) :: T(3,3)
      real(8), intent(inout) :: EE(3,3,3,3)
      real(8) ::  p ,fac,nn
       
      p = -tr(T)/3.0d0
      nn =  mat%xismallE - mat%nsmallE
      fac = (1.0d0 - exp(-p/mat%prefsmallE))**(nn)
      if (fac < 0.005) return              
      !pref = 1.0d0
      !xi = 2.0d0
      !n = mat%n
      !pros = -tr(Tb)/3.0d0            
      EE = EE  *fac
      end subroutine smallpressurestiffness
      !-----------------------------------------------------------------  
      
      
      



      end subroutine constitutiveRates
      !-----------------------------------------------------------------
      
      
      
      end module wo2_sch
      !-----------------------------------------------------------------