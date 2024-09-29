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

!=================================================================
      module umat_tools
      implicit none
      integer :: rnoel,rnpt,rkstep,rkinc,rndi,rnshr
      real(8),allocatable ::rstress(:),rstatev(:),rdstran(:),rprops(:)
      real(8) ::  rtime(2),rdtime,rcoords(3)  
      character*80 :: rcmname
      
      
      !=================================================================
      contains
      !=================================================================
      
      !-----------------------------------------------------------------
      subroutine set_report(stress,statev,dstran,time,dtime, cmname,
     &    ndi,nshr,ntens, nstatv,props,nprops,coords, noel,npt, 
     &    kstep,kinc)
      integer :: ntens,nstatv,nprops,noel,npt, kstep,kinc,ndi,nshr
      real(8) ::stress(ntens),statev(nstatv),dstran(ntens),
     &   time(2),dtime,props(nprops),coords(3) 
      character*80 :: cmname
      
      rcmname = cmname
      rnoel=noel; rnpt=npt; rkstep=kstep; rkinc=kinc; 
      rndi=ndi;rnshr=nshr
      rtime = time; rdtime = dtime
      if (.not. allocated(rstress)) allocate(rstress(ntens))
      if (.not. allocated(rstatev)) allocate(rstatev(nstatv))
      if (.not. allocated(rdstran)) allocate(rdstran(ntens))
      if (.not. allocated(rprops)) allocate(rprops(nprops))
      rstress(:) = stress(:)
	  rstatev(:) = statev(:); rdstran(:) = dstran(:)
      rprops(:) = props(:); rcoords(:) = coords(:) 
      end subroutine set_report
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      subroutine report(msg)      
      implicit none
      character(len=*),intent(in) :: msg ! when passing strings as arguments of a subroutine, the string should be declared with (len=*)            
      integer,dimension(8) :: dat ! arguments for date_and_time
      integer :: i,j,unit      
      unit=6
      write(*,*)  'msg=',msg 
      write(*,*)  'Details about this error in unit ',unit 
      !call reportFil(msg)      
      write(unit,*) '********************************'
      write(unit,*) 'MATERIAL NAME = ', trim(rcmname)       
      write(unit,*) msg      
      write(unit,*) '********************************'
      call date_and_time(values=dat)      
      write(unit,'(A,2(I2.2,A),I4.4,6X,A,I2.2,A,I2.2,A,I2.2,A,I3.3)')
     &   ' DATE: ',dat(3),'.', dat(2),'.', dat(1), 'TIME: ',
     &    dat(5),':',dat(6),':',dat(7),'.',dat(8)
      write(unit,*) 'UMAT was called with the following arguments:'
      write(unit,'(4(A,I7))') ' kstep=',rkstep,' kinc=',rkinc       
      write(unit,'(2(A,I8))')' ELEMENT No.=',rnoel,' GAUSS POINT=',rnpt
      write(unit,*) 'COORDS='
      write(unit,'(3(es15.7e3,2X) )')  rcoords
      write(unit,'(2(A,g13.5,2X))') ' STEP TIME=',rtime(1),
     &    'TOTAL TIME=',rtime(2)
      write(unit,'(2(A,g13.5,2X))') ' TIME INCREMENT=',rdtime 
      write(unit,*) 'MATERIAL PROPERTIES:'
      do i=1, size(rprops)
        write(unit,'( A8,i2,A2,2X,g13.4e3)') '  props(',i,')=',
     &              rprops(i)
      end do
      write(unit,*) 'TOTAL STRESS stress='
      write(unit,'(6(es15.7e3,2X) )')  rstress
      write(unit,*) 'STRAIN INCREMENT dstran='
      write(unit,'(6(es15.7e3,2X) )')  rdstran
      write(unit,*) 'STATE VARIABLES statev='
      write(unit,'(5(es15.7e3,2X) )')  rstatev
      
      if('ERR' .ne. msg(1:3))  return 
      call xit   
      end subroutine report
	  !----------------------------------------------------------------- 

           
      !-----------------------------------------------------------------  
      subroutine elastic(stress,ddsdde,dstran,ntens,E,nu)
      use niemunis_tools_nano
      implicit none
      integer,intent(in) :: ntens      
      real(8),intent(in) :: dstran(ntens),E,nu
      real(8),intent(inout) :: stress(ntens)
      real(8),intent(out) :: ddsdde(ntens,ntens)
      real(8) :: Tb(3,3),cE(3,3,3,3),depsb(3,3)            
      
      Tb = map2T(stress,ntens)
      depsb = map2D(dstran,ntens)       
      cE = IsoElastic(E=E,nu=nu)
      Tb = Tb + (cE .xx. depsb)
      stress = map2stress(Tb,ntens)
      ddsdde = map2ddsdde(cE,ntens) 
      end subroutine elastic
      !-----------------------------------------------------------------
      
       !--------------------------------------------------------------------- 
       subroutine checkRange(x,xmin,xmax,xname)
       implicit none
       real(8),intent(in) :: x,xmin,xmax
       character(len=*),intent(in) ::  xname
       character(len=40) :: msg      
       if( x< xmin .or. x> xmax) then
         msg = "ERROR props2mat: " // trim(xname) // ' out of range'
         call report(msg)
       end if
       end subroutine checkRange
       !--------------------------------------------------------------------- 
      
!      !------------------------------------------------------------------------
!      Pure Function to_upper (str) Result (string)
!      !   ==============================
!      !   Changes a string to upper case
!      !   ==============================
!      Implicit None
!      Character(*), Intent(In) :: str
!      Character(LEN(str))      :: string
!      Integer :: ic, i
!      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!      Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
!      !   Capitalize each letter if it is lowecase
!      string = str
!      do i = 1, LEN_TRIM(str)
!        ic = INDEX(low, str(i:i))
!        if (ic > 0) string(i:i) = cap(ic:ic)
!      end do
!      End Function to_upper
      !------------------------------------------------------------------------   

      !------------------------------------------------------------------------   
      subroutine to_upper(str,string)      
      Implicit None
      Character(*), Intent(In) :: str
      Character(LEN(str)),intent(out)      :: string
      Integer :: ic, i
      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
      string = str
      do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
      end do
      end subroutine to_upper
      !------------------------------------------------------------------------   
      
      end module umat_tools
      !=================================================================