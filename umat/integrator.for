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
	  
	  module integrator
      use umat_tools      
      use wo2_sch  
      implicit none
            
      
      contains          
       
      
      
      !-----------------------------------------------------------------
      subroutine rungeKutta(y,D,props,ntens,nQ,nprops,dtime,TOLT,ATOL,
     &      TOLTred,nMaxSubi,nMaxIter,asv,nasv,pnewdt,hTot)
      use niemunis_tools_nano, only: norm
      implicit none
      integer, intent(in) :: ntens,nQ,nprops,nMaxSubi,nMaxIter,nasv
      real(8), intent(in) :: D(3,3), props(nprops), dtime, TOLT,ATOL,
     &      TOLTred
      real(8), intent(inout) :: y(ntens+nQ),
     &      pnewdt, hTot
      real(8), intent(out) :: asv(nasv)
      real(8),dimension(ntens+nQ) :: dydt0,v,ym,
     &      dydtm,w,ynew,dydtnew
      real(8) :: h, fhnew,errorT,TOL, dEpsmax
      integer :: errCode, nSubi,isubi,iiter      
      
      dEpsmax = 1.0d-5    ! used for the first increment only
      nSubi= max(int(norm(D)*dtime/dEpsmax),1)
      h = dtime/nSubi   ! start with an small time increment for better integration of interg. strain  
      !h = dtime*0.1d0 ! this may have a large impact on the convergence in Abaqus      
      hTot = 0
      fhnew = 1.0d0
      if(y(1)>-1.0d0) then
      continue
      end if
      ! first call to the evolution equations       
      call getdydt(y,D,h,props,ntens,nQ,nprops,dydt0,asv,nasv,errCode)      
!      if(dydt0(1)==0) dydt0(1) = -1   ! if T is on the bounding, dTdt =0 and stress cannot evolve    
      if(errCode == 1) call report('ERROR rungeKutta: bad y')
!     &     // 'equation cannot be evaluated for the current state y')
      do isubi=1,nmaxsubi
        call getTOL(TOLT,TOLTred,errCode,TOL)         
        do iiter=1,nmaxiter          
          if(h<1.0d-10 .and. h < dtime-hTot) then
            write(*,*)! 'rnoel',rd%rnoel, 'rnpt',rd%rnpt,
     &      'h is too small wo2'
!            pnewdt = 0.5d0 ! abaqus will reattempt the increment with an smaller time increment
            return
          end if
          ym = y + dydt0* h/2.0d0
          v =  y + dydt0* h 
          call getdydt(ym,D,h,props,ntens,nQ,nprops,dydtm,asv,nasv,
     &                  errCode) ! evaluate the constitutive equations for the mid. of increment
          if(errCode == 1) then  ! problems evaluating the constituve equation
            h = h/2.0d0
            cycle ! try again with an smaller h
          end if
          w = ym + dydtm*h/2.0d0
          call getErrorT(v(1:ntens),w(1:ntens),y(1:ntens),ntens,
     &            ATOL,errorT)
          if(errorT == 0) fhnew = 5.0d0
          if(errorT > 0) fhnew = min(5.0d0,0.9d0*sqrt(TOL/errorT)) 
!          call saveRungeKuttaInfo(isubi,iiter,hTot,h,dtime,y(1:ntens),
!     &      ym(1:ntens),v(1:ntens),w(1:ntens),
!     &      dydt0(1:ntens),dydtm(1:ntens),errorT,TOL,ntens) ! for debug only
          if(errorT < TOL) exit          
          h = max(0.1d0, fhnew) * h
        end do ! end iterations
        if(iiter >= nmaxiter) then            
           if(dydt0(1)==-1) dydt0(1) = 0  ! Trying to go outside bounding: continue on the bounding           
           h = h/2.0d0
           cycle ! try again with an smaller h           
        end if
        ynew = 2.0d0*w - v  ! new solution
        call getdydt(ynew,D,h,props,ntens,nQ,nprops,dydtnew,asv,nasv,
     &               errCode) ! check if the new solution can be evaluated by the const. equation
        if(errCode == 1) then
          h = h/2.0d0
          cycle ! try again with an smaller h
        endif
        ! accept the solution
        y = ynew
        dydt0 = dydtnew
        hTot = hTot + h
        asv(1) = isubi
        !qav0 = qav
        if(hTot >= dtime) exit
        h = max(0.2d0, fhnew) * h
        h = min(h, dtime-hTot)       
      end do ! end subincrements
      if(isubi >= nmaxsubi) then  
  !      write(*,*) ! 'rnoel',rd%rnoel, 'rnpt',rd%rnpt,
   !  &   'too many subincrements required'
   !        pnewdt = 0.5d0 ! abaqus will reattempt the increment with an smaller time increment
   !     call report('WARNING rungeKutta: '
   !  &      //      ' too many subincrements are required')
        return
      end if
      end subroutine rungeKutta
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine getdydt(y,D,dt,props,ntens,nQ,nprops,dydt,asv,nasv,
     &      errCode)      
      use niemunis_tools_nano, only: map2stress,map2T,norm
      implicit none
      integer,intent(in) :: nQ,nprops,ntens,nasv
      real(8),intent(in) :: y(ntens+nQ),D(3,3),
     &  props(nprops),dt
      real(8),intent(out) :: dydt(ntens+nQ),
     &     asv(nasv)
      integer,intent(out) :: errCode
      real(8) :: T(3,3),dTdt(3,3),Q(nQ),dQdt(nQ),B(ntens,ntens),
     &  Tvar(3,3),Dvar(3,3),Tvec(ntens),Qvar(nQ),G(nQ,ntens),
     &  dTdtvar(3,3),dQdtvar(nQ),dBdt(ntens,ntens),AuxT(3,3),theta,
     &  normD,dGdt(nQ,ntens),minD
      integer, parameter,dimension(1:6) :: i6=(/ 1,2,3,1,1,2/),        ! table of index for transformations Tensor-Matrix
     &                                     j6=(/ 1,2,3,2,3,3/)
      integer :: m,iB,iG,iQ,i,j,errCodeJ

      errCode = 0  ! everything ok
      dydt(:) = 0
      normD = norm(D)
      theta = max(normD,1.d0)*1.0d-7
      call getMinNonZeroVal(D,3,minD)
      if (minD>0 .and. theta > 1.0d-3*minD)  theta = 1.0d-3*minD     
      
      T = map2T(y(1:ntens),ntens)      
      Q = y(ntens+1:ntens+nQ)      
      call constitutiveRates(T,Q,D,dt,props,ntens,nQ,nprops,dTdt,dQdt,
     &      asv,nasv,errCode,0)        
      if(errCode == 1) return ! error in evaluating the constitutive equation (T, or Q inadmissible or props not ok)
      dydt(1:ntens) = map2stress(dTdt,ntens)  ! dTdt to y 
      dydt(ntens+1:ntens+nQ) = dQdt  ! dQdt to y	

      end subroutine getdydt
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine get_ddsdde0(ddsdde,y,D,props,ntens,nQ,nprops,dtime,
     &                      asv,nasv)
      use niemunis_tools_nano, only: map2stress,map2T,norm
      implicit none
      integer, intent(in) :: ntens,nasv,nQ,nprops
      real(8), intent(in) :: y(ntens+nQ),D(3,3),props(nprops),dtime 
      real(8), intent(out) :: asv(nasv) 
      real(8), intent(out) :: ddsdde(ntens,ntens)
      real(8) :: T(3,3),dTdt(3,3),dTdt1(3,3),dTdt2(3,3), D1(3,3),
     &      D2(3,3), Q(nQ),dQdt(nQ),  normD, theta
      integer :: errCod 
      integer, parameter,dimension(1:6) :: i6=(/ 1,2,3,1,1,2/),        ! table of index for transformations Tensor-Matrix
     &                                     j6=(/ 1,2,3,2,3,3/)
      integer :: m,i,j
      
      normD = norm(D)
      theta = max(normD,1.d0)*1.0d-7
      
      T = map2T(y(1:ntens),ntens)
      Q = y(ntens+1:ntens+nQ)
      
      do m=1,ntens
        i = i6(m)
        j = j6(m)
        D1 = 0
        D2 = 0
        if(m<=3) then        
          D1(i,j) = D1(i,j) - theta
          D2(i,j) = D2(i,j) + theta
        end if
        if(m>3) then
          D1(i,j) = D1(i,j) - theta/2.0d0
          D1(j,i) = D1(j,i) - theta/2.0d0
          D2(i,j) = D2(i,j) + theta/2.0d0
          D2(j,i) = D2(j,i) + theta/2.0d0
        end if
        call constitutiveRates(T,Q,D1,dtime,props,ntens,nQ,nprops,
     &        dTdt1,dQdt,asv,nasv,errCod,1)
        call constitutiveRates(T,Q,D2,dtime,props,ntens,nQ,nprops,
     &         dTdt2,dQdt,  asv,nasv,errCod,1)
        dTdt = (dTdt2 - dTdt1)/(2.0d0*theta)
        ddsdde(:,m) = map2stress(dTdt,ntens )
      end do
      end subroutine get_ddsdde0
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine get_ddsdde(ddsdde,y,D,props,ntens,nQ,nprops,dtime,
     &                      asv,nasv)
      use niemunis_tools_nano, only: map2stress,map2T,norm
      implicit none
      integer, intent(in) :: ntens,nasv,nQ,nprops
      real(8), intent(in) :: y(ntens+nQ),D(3,3),props(nprops),dtime 
      real(8), intent(out) :: asv(nasv) 
      real(8), intent(out) :: ddsdde(ntens,ntens)
      real(8) :: T(3,3),dTdt(3,3),dTdt1(3,3),dTdt2(3,3), D1(3,3),
     &      D2(3,3), Q(nQ),dQdt(nQ),  normD, theta,dTdtdtheta(3,3)
      integer :: errCod 
      integer, parameter,dimension(1:6) :: i6=(/ 1,2,3,1,1,2/),        ! table of index for transformations Tensor-Matrix
     &                                     j6=(/ 1,2,3,2,3,3/)
      integer :: m,i,j
      
      normD = norm(D)
      theta = max(normD,1.d0)*1.0d-7
      
      T = map2T(y(1:ntens),ntens)
      Q = y(ntens+1:ntens+nQ)
      
      call constitutiveRates(T,Q,D,dtime,props,ntens,nQ,nprops,
     &        dTdt,dQdt,asv,nasv,errCod,1)
      
      do m=1,ntens
        i = i6(m)
        j = j6(m)
        D2 = D
        if(m<=3) then                  
          D2(i,j) = D(i,j) + theta
        end if
        if(m>3) then           
          D2(i,j) = D(i,j) + theta/2.0d0
          D2(j,i) = D(j,i) + theta/2.0d0
        end if
        call constitutiveRates(T,Q,D2,dtime,props,ntens,nQ,nprops,
     &         dTdt2,dQdt,  asv,nasv,errCod,1)
        dTdtdtheta = (dTdt2 - dTdt)/(theta)
        ddsdde(:,m) =  map2stress(dTdtdtheta,ntens )
      end do
      end subroutine get_ddsdde
      !-----------------------------------------------------------------
      
      
      
      
      !-----------------------------------------------------------------
      subroutine getJacobiD0(stress,statev,ddsdde,props,dtime,theta,
     &      ntens,nstatev,nQ,nprops)
      use niemunis_tools_nano, only: map2stress,map2T
      implicit none
      integer, intent(in) :: ntens,nstatev,nQ,nprops
      real(8), intent(in) :: stress(ntens),statev(nstatev), 
     &      props(nprops),dtime,theta
      real(8), intent(out) :: ddsdde(ntens,ntens)
      real(8) :: T(3,3),dTdt(3,3),dTdt1(3,3),dTdt2(3,3), D1(3,3),
     &      D2(3,3), Q(nQ),dQdt(nQ), dsdt(ntens), asv(nstatev-nQ)
      integer :: errCod,nasv
      integer, parameter,dimension(1:6) :: i6=(/ 1,2,3,1,1,2/),        ! table of index for transformations Tensor-Matrix
     &                                     j6=(/ 1,2,3,2,3,3/)
      integer :: m,i,j
      
      T = map2T(stress,ntens)
      Q = statev(1:nQ)
      nasv = size(asv)
      if (nasv > 0) asv(:) = statev(nQ+1:nQ + nasv)

      do m=1,ntens
        i = i6(m)
        j = j6(m)
        D1 = 0
        D2 = 0
        if(m<=3) then        
          D1(i,j) = D1(i,j) - theta
          D2(i,j) = D2(i,j) + theta
        end if
        if(m>3) then
          D1(i,j) = D1(i,j) - theta/2.0d0
          D1(j,i) = D1(j,i) - theta/2.0d0
          D2(i,j) = D2(i,j) + theta/2.0d0
          D2(j,i) = D2(j,i) + theta/2.0d0
        end if
        call constitutiveRates(T,Q,D1,dtime,props,ntens,nQ,nprops,
     &        dTdt1,dQdt,asv,nasv,errCod,1)
        call constitutiveRates(T,Q,D2,dtime,props,ntens,nQ,nprops,
     &         dTdt2,dQdt, asv,nasv,errCod,1)
        dTdt = (dTdt2 - dTdt1)/(2.0d0*theta)
        ddsdde(:,m) = map2stress(dTdt,ntens )
      end do
      end subroutine getJacobiD0
      !-----------------------------------------------------------------
      
      

      
      !-----------------------------------------------------------------
      subroutine getErrorT(v,w,y,ntens,ATOL,errorT)
      implicit none
      integer, intent(in) :: ntens
      real(8), intent(in) :: v(ntens),w(ntens),y(ntens),ATOL
      real(8), intent(out) :: errorT
      real(8) :: aux(ntens),s,ymax,rmax,rmin,wy(ntens)
      integer :: i
      
      ! original definition of error after Fellin.
      errorT = 0
      do i=1,ntens
        s = max(abs(y(i)),abs(w(i))) + ATOL    ! add 1 to avoid small numbers in the divisor??
        errorT = errorT + ((w(i) - v(i))/s)**2
      end do
      errorT = sqrt(errorT)
      end subroutine getErrorT 
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      subroutine getTOL(TOLT,TOLTred,errCode,TOL) ! returns the tolerance in stress
      implicit none                            ! if errCode=2 (i.e. for intergranular strain)
      integer, intent(in) :: errCode           ! then the reduced tolerance TOLTred will be used
      real(8), intent(in) :: TOLT, TOLTred
      real(8), intent(out) :: TOL

      if(errCode == 0) TOL = TOLT    ! normal tolerance
      if(errCode == 2) TOL = TOLTred ! reduced tolerance tolerance (intergranular strain)     
      end subroutine getTOL
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine forwardEuler(y,D,props,ntens,nQ,nprops,dtime,dEpsmax,
     &   asv,nasv)
      use niemunis_tools_nano
      implicit none
      integer, intent(in) :: ntens,nQ,nprops,nasv
      real(8), intent(in) :: D(3,3), props(nprops), dtime, dEpsmax
      real(8), intent(inout) :: y( ntens+nQ)
      real(8), intent(out) :: asv(nasv)
      real(8),dimension( ntens+nQ) :: dydt
      real(8) :: h 
      integer,parameter :: nMaxSubi = 10000
      integer ::  nSubi,errCode,isubi
      
      nSubi= max(int(norm(D)*dtime/dEpsmax),1)
      h = dtime/nSubi
      do isubi=1,nSubi
        call getdydt(y,D,h,props,ntens,nQ,nprops,dydt,asv,nasv,errCode)
        if(errCode == 1) then
          call report('ERROR in umat forwardEuler: bad state ')
!     &      //'equation could not be evaluated. '
!     &      //'T or state out of range') 
          call xit
        end if
        y = y + dydt*h
        asv(1) = isubi          
      end do ! end subincrements
      end subroutine forwardEuler
      !-----------------------------------------------------------------
      
      
      !-----------------------------------------------------------------
      subroutine getMinNonZeroVal(A,n,val)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: A(n,n)
      real(8), intent(out) :: val
      real(8) :: x
      integer :: i,j
      val = 1.0d99
      do i=1,n
        do j=1,n
          x = abs(A(i,j))
          if ( x > 0.0d0 .and. x < val) val = x
        end do
      end do      
      end subroutine getMinNonZeroVal
      !-----------------------------------------------------------------
      
      end module integrator