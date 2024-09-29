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
      subroutine UMAT(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      use niemunis_tools_nano, only: map2D
      use wo2_sch
      use shc_frac, only: get_Eampl2,get_dNdt
      use integrator
      use umat_tools
      implicit none     
      
      character*80 cmname, mname 
      character*256  reportfName 
      integer ndi,nshr,ntens,nstatv,nprops,npt,layer,kspt,kstep,kinc,
     &        noel
      double precision ::  stress(ntens),statev(nstatv),
     &       ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &       stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     &       props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &       dtime,sse,spd,scd,drpldt,temp,dtemp,rpl,pnewdt,celent
      ! End of variable declaration for umat
       
      integer :: nQ   ! number of state variables: epor + pwater + 6 of intergranular strain
      integer :: nasv      ! =nstatv - nQ : number of additional state variables (used for output purposes only)
      integer :: nMaxSubi  ! max. number of Subincrements in RungeKutta
      integer :: nMaxIter  ! max. number of iterations in RungeKutta
      real(8) :: T(3,3),D(3,3),dTdt(3,3),B(ntens,ntens), 
     &    TOLT,ATOL,theta,normD,depsMax,TOLTred, E,nu,
     &    Baux(ntens,ntens),dt,hTot,hTry,ddsdde0(ntens,ntens)
      real(8), allocatable :: Q(:),dQdt(:),y(:)
      real(8), allocatable :: asv(:)  ! additional state variables used for output purposes only (different from Q)
      integer :: integrMethod ! 1-> RungeKutta, 2-> ForwardEuler
      logical :: isElas
      real(8) :: Eampl,dNdt ! for HCA
      
      call set_report(stress,statev,dstran,time,dtime, cmname,         ! copy all variables to be printed in report
     &    ndi,nshr,ntens, nstatv,props,nprops,coords, noel,npt, 
     &    kstep,kinc) 
      
      if(any(isnan(statev)).or. any(isnan(stress)).or.                              ! \com check if any  NANs in the input
     &   any(isnan(dstran)) ) then
           call report('ERROR: there is nan in the input')
      end if
  
      !----------------------------------------------------------------
      ! PARAMETERS FOR INTEGRATION METHOD
      !----------------------------------------------------------------
      integrMethod = 1  ! integration Method: 1-> RungeKutta, 2-> ForwardEuler
      TOLT = 1.0d-3     ! tolerance in the stress
      ATOL = 1.0d-3     ! absolute tolerance (to avoid null division)
      TOLTred = 1.0d-4 ! reduced stress tolerance (i.e. for a reversal in intergranular strain or plasticity)       
      nMaxSubi = 1000   ! max. number of Subincrements in RungeKutta
      nMaxIter = 20     ! max. number of iterations/attempts in RungeKutta
      depsMax = 1.0d-5  ! maximal size of strain increment for forward euler  
      !----------------------------------------------------------------      

	!----------------------------------------------------------------
      ! JACOBIAN      
	!ctrl%computeJacobian = .true. ! true for abaqus/standar, false for explicit
	!----------------------------------------------------------------             
       
      !----------------------------------------------------------------      
      isElas = .false.   ! True to compute with isotropic elasticity
      !if (kstep==1)  isElas = .true.
      if(isElas) then
         E =  1.0d6; nu = 0.3d0;
         call elastic(stress,ddsdde,dstran,ntens,E,nu)
         return
      end if      
      !----------------------------------------------------------------      
      
      !call get_Eampl(Eampl,KSTEP,KINC,TIME,NOEL,NPT,COORDS)
      !call get_dNdt(dNdt,KSTEP,KINC,TIME,NOEL,NPT,COORDS)
      !statev(28:29) = [Eampl,dNdt]
      
      
      call get_nQ_nasv(nQ,nasv,ntens,nstatv)
      allocate(y(ntens+nQ),Q(nQ),asv(nasv))
            
      call statev2Q(statev,nstatv,Q,nQ,ntens) ! state variables with evolution equations
      !call statev2asv(statev,nstatv,asv,nasv)
      
      
      y(1:ntens) = stress(:)
      y(ntens+1:ntens+nQ) = Q(:)          
      
      if (kstep==2) then
        !call report('ERR second step')
      end if 

      if(dtime.le.1.0d-40)  D  = map2D(dstran,ntens)                    !   handle the case $\Db=0$,dtime=0
      if(dtime.gt.1.0d-40)  D  = map2D(dstran,ntens)/dtime 
      
      !if((dtime<1.0d-40).and.(dot_product(dstran,dstran)<1.0d-40)) then  ! in Abaqus 6.14 dtime =1 for first increment of first step
   !   if ((kstep*kinc .le. 1) .or. 
   !  &      (dot_product(dstran,dstran)<1.0d-40) .or.
   !  &      (dtime<1.0d-40)     ) then
      !if  ((kstep*kinc .le. 1) ) then
       !integrMethod = 2 ! explicit Euler for the first increment?? (geostatic)
       !!call get_ddsdde0(ddsdde,y,D,props,ntens,nQ,nprops,dtime,asv,nasv)       
       !!return      ! only jacobian required, no stress update
      !end if               
      select case(integrMethod)
        case(1)          
            call rungeKutta(y,D,props,ntens,nQ,nprops,dtime,TOLT,ATOL,
     &              TOLTred,nMaxSubi,nMaxIter,asv,nasv,pnewdt,hTry)                                    
        case(2)
            call forwardEuler(y,D,props,ntens,nQ,nprops,dtime,depsMax,
     &                      asv,nasv)            
      end select      
      call get_ddsdde(ddsdde,y,D,props,ntens,nQ,nprops,dtime,asv,nasv)  
      
      stress(:) = y(1:ntens)   
      Q(:) = y(ntens+1:ntens+nQ) 
      call Q2statev(Q,nQ,statev,nstatv,ntens)
      call asv2statev(asv,nasv,statev,nstatv)      
      !is = map2D(statev(3:8),6)  ! intergranular strain
      !call rigidrot3(is, dfgrd0,dfgrd1)  
      !statev(3:8) = map2stran(is,ntens) 
      
      if(any(isnan(statev)).or.   any(isnan(stress)).or.                            ! \com check if any  NANs in result
     &   any(isnan(ddsdde))     )then        
        write(*,*) 'stress'
        write(*,*) stress
        write(*,*) 'statev'
        write(*,*) statev
        write(*,*) 'ddsdde'
        write(*,*) ddsdde      
        call report('ERROR umat: there is nan in the output')
      end if      
      
      !--------------------
      !contains
      !--------------------
      
      
      
      
      end subroutine UMAT 
      !-----------------------------------------------------------------
      
      
      
      