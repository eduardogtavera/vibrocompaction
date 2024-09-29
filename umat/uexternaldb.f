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
	  
	  
	  
	  !------------------------------------------------------------------------
      ! This file is used to set up the sequence and the velocity in which
      ! the perforations will take place. At the beginning of the compuation
      ! a text file with this information will be read and their values saved 
      ! in a global module that can be accesed from initial conditions and umat.
      ! 08.03.2023 C.Grandas
      !------------------------------------------------------------------------
      module globalVariables
      implicit none            
      character(256) ::  file_initial,fcheck
      integer :: ufcheck
      logical :: printPileSeq
      integer :: intv(500)
      real(8) :: REALV(500)
      character(8) :: CHARV(500)
      integer :: nbhole, nInstStep
      integer, dimension(:), allocatable :: inst_sequence
	  !integer, allocatable(:) :: bhole_label
      real(8), dimension(:,:), allocatable :: bhole_top, bhole_bot
      real(8), dimension(:), allocatable :: freq, inst_velocity      
      !file_initial="C:/Data/projects/InputBaugrube3D.txt"     
      end module globalVariables
      !------------------------------------------------------------------------

      !----------------------------------------------------------------
      SUBROUTINE readBoreInfo( )    
      ! read data about installation sequence.  
      use globalVariables
      implicit none 
      integer :: ufile, klabel, i, kstep, kbhole
      character(500) :: line
	  real(8) :: xtop, ytop
      
      ufile = 54  
      open (unit=ufile, file=file_initial,ACTION='READ', err=333) 
      read(ufile,*) line ! *nbhole
      read(ufile,*) nbhole
      allocate(bhole_top(nbhole*2,3),bhole_bot(nbhole*2,3) ) ! alocate large arrays
      allocate(freq(nbhole),inst_velocity(nbhole) )
	  !allocate(bhole_label(nbhole))

      read(ufile,*) line ! *xtop, ytop, ztop, xbot, ybot, zbot, freq, inst_velocity
      do i=1,nbhole
        read(ufile,*) klabel, xtop, ytop
		bhole_top(klabel,1:2) = [xtop, ytop] !, bhole_bot(i,1:3), freq(i), 
!     &        inst_velocity(i)
      end do
      
      read(ufile,*) line ! *nInstSteps
      read(ufile,*) nInstStep
      allocate(inst_sequence(2*nInstStep)) ! alocate a large array
      inst_sequence(:) = 0
      read(ufile,*) line !  *step, bholes 
      do i=1,nInstStep
        read(ufile,*)  kstep, kbhole
		inst_sequence(kstep) = kbhole
      end do
      close (ufile)
      return

  333 CALL STDB_ABQERR(-3,'Error in uexternaldb:  '
     &   //   'I cannot open the file '//trim(file_initial),
     &    intv,REALV,CHARV) ! to stop the calculation!	
      
      end subroutine readBoreInfo
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
      
      ! user coding to set up the FORTRAN environment, open files, close files, 
      ! calculate user-defined model-independent history information,
      ! write history information to external files,
      ! recover history information during restart analyses, etc.
      ! do not include calls to utility routine XIT
            
      !----------------------------------------------------------------------     
      !  The use of ABA_PARAM.INC eliminates the need to have different
      !  versions of the code for single and double precision.  
      !  ABA_PARAM.INC defines an appropriate IMPLICIT REAL statement and 
      !  sets the value of NPRECD to 1 or 2, depending on whether the 
      !  machine uses single or double precision.              
      !----------------------------------------------------------------------     
      !
      !  Variables:  (all variables are passed in for information. No 
      !               variables need to be defined.)
      !  LOP
      !
      !    LOP=0 indicates that the subroutine is being called at the start of 
      !          the analysis.
      !    LOP=1 indicates that the subroutine is being called at the start of 
      !          the current analysis increment. 
      !    LOP=2 indicates that the subroutine is being called at the end of 
      !          the current analysis increment. When LOP=2, all information 
      !          that you need to restart the analysis should be written to 
      !          external files.
      !    LOP=3 indicates that the subroutine is being called at the end of 
      !          the analysis.
      !    LOP=4 indicates that the subroutine is being called at the beginning 
      !          of a restart analysis. When LOP=4, all necessary external 
      !          files should be opened and properly positioned and all 
      !          information required for the restart should be read from 
      !          the external files.
      !  
      !  LRESTART
      !
      !    LRESTART=0 indicates that an analysis restart file is not being 
      !               written for this increment.
      !    LRESTART=1 indicates that an analysis restart file is being 
      !               written for this increment.
      !    LRESTART=2 indicates that an analysis restart file is being 
      !               written for this increment and that only one increment 
      !               is being retained per step so that the current increment 
      !               overwrites the previous increment in the restart file.
      !
      !  TIME =    Current (step,total) time.!    
      !  DTIME =   Time increment    
      !  KSTEP =   Current step number. When LOP=4, KSTEP gives the restart 
      !                 step number.
      !  KINC =    Current increment number. When LOP=4, KINC gives the 
      !              restart increment number.
      !
      !  Other variables
      !  FNAME =   Root file (without extension) name of input file(s).
      !  OUTDIR =  Path of the current analysis directory

      
      use globalVariables
      !INCLUDE 'ABA_PARAM.INC'  
      IMPLICIT NONE      
      REAL(8) :: TIME(2),DTIME
      INTEGER :: LOP,LRESTART,KSTEP,KINC
      
      !DIMENSION TIME(2)
         
      character (256) ::  outdir,jobname
      integer :: lenoutdir,lenjobname
      integer,dimension(8) :: dat ! arguments for date_and_time
      character(1) :: separator
	  integer :: kpile
	  real(8) :: xtop,ytop
      
      write(*,*) ''	
      write(*,*) 'STARTING UEXTERNALDB kstep=',kstep,'kinc=',kinc
      write(*,*) 'lop',lop	
      
      CALL STDB_ABQERR(1,'uexternaldb called at '
     &         // 'step=%I, increment=%I with lop=%I and lrestart=%I.',
     &           [kstep,kinc,lop,lrestart],REALV,CHARV)
      
      select case(lop)
        case(0) 
        ! LOP=0 indicates that the subroutine is being called at 
        !       the start of the analysis. 
        call getoutdir( outdir, lenoutdir )

        ! choose '\' or '/' as folder separation depending on OS
        !DEC$ IF DEFINED(_WIN32)
            separator = '\'
        !DEC$ ELSEIF DEFINED(__linux)
            separator = '/'
        !DEC$ ELSE
            separator = '/'
        !DEC$ ENDIF        
        
        file_initial = trim(outdir) //separator//'piles_info.txt'  
        !file_initial = trim(outdir) // '\' // 'InputBaugrube3D.txt' ! windows
        !file_initial = trim(outdir) // '/' // 'InputBaugrube3D.txt' ! linux
        call readBoreInfo( ) 
		
        ! check pile installation sequence
		fcheck = trim(outdir) //separator//'checkInstalation.txt'  
        ufcheck = 44        
        open (unit=ufcheck, file=fcheck,err=43)  
        printPileSeq = .true.		
        
        CALL STDB_ABQERR(1,'uexternaldb called at the start of the '
     &    // 'analysis. Layer and geometry information '
     &    // 'obtained '
     &    // 'from: '//trim(file_initial) // '. ',
     &    intv,REALV,CHARV)
        CALL STDB_ABQERR(1,'Installation sequence written in: '
     &    // trim(fcheck) // '. ',
     &    intv,REALV,CHARV)
        return                                                      
  43    CALL STDB_ABQERR(-3,'Error in uexternaldb:  '
     &   //   'I cannot open the file '//trim(fcheck),
     &    intv,REALV,CHARV) ! to stop the calculation! 	 
     

      case(1) 
        !    LOP=1 indicates that the subroutine is being called at the start of 
        !          the current analysis increment.
        CALL STDB_ABQERR(1,'uexternaldb called at the start of '
     &         // 'the current increment.',intv,REALV,CHARV)  
	 
		kpile = inst_sequence(kstep)
		if (kpile < 1) return
		xtop = bhole_top(kpile,1)
		ytop = bhole_top(kpile,2)
		write(*,*) 'kstep=', kstep
		write(*,*) 'kpile=', kpile, 'xtop=', xtop, 'ytop=',ytop
		
	 
        !call timeStamp(kstep,kinc,'INC START')
        return
		
      case(2) 
        !    LOP=2 indicates that the subroutine is being called at the end of 
        !          the current analysis increment. When LOP=2, all information 
        !          that you need to restart the analysis should be written to 
        !          external files.
        printPileSeq = .false.
		
        CALL STDB_ABQERR(1,'uexternaldb called at the end of '
     &       // 'the current increment.', intv,REALV,CHARV) 
        !call timeStamp(kstep,kinc,'INC END')  
        return
                        
      case(3) 
        !    LOP=3 indicates that the subroutine is being called at the end of 
        !          the analysis.		
         CALL STDB_ABQERR(1,'uexternaldb called at the end of '
     &         // 'the analysis.', intv,REALV,CHARV)    
         return      

      case(4)
        ! LOP=4 indicates that the subroutine is being called at the beginning 
        !          of a restart analysis. When LOP=4, all necessary external 
        !          files should be opened and properly positioned and all 
        !          information required for the restart should be read from 
        !          the external files.          

        call getoutdir( outdir, lenoutdir )         
        !DEC$ IF DEFINED(_WIN32)
            separator = '\'
        !DEC$ ELSEIF DEFINED(__linux)
            separator = '/'
        !DEC$ ELSE
            separator = '/'
        !DEC$ ENDIF                
        
        !fcheck = trim(outdir) //separator//'stepPileCoords.txt' ! check pile installation sequence 
        !ufcheck = 44        
        !open (unit=ufcheck, file=fcheck,err=44)
        printPileSeq = .true.

        
        CALL STDB_ABQERR(1,'uexternaldb called at the beginning of '
     &         // 'a restart analysis.',      intv,REALV,CHARV)
        return
  44    CALL STDB_ABQERR(-3,'Error in uexternaldb:  '
     &   //   'I cannot open the file '//trim(fcheck),
     &    intv,REALV,CHARV) ! to stop the calculation!        
      end select    

      RETURN            
      
      contains        
      
      !----------------------------------------------------------------
      subroutine timeStamp(kstep,kinc,label)
      implicit none
      character(*),intent(in):: label
      integer, intent(in):: kstep,kinc
      integer,dimension(8) :: dat ! arguments for date_and_time
      
      call date_and_time(values=dat)
      write(191,'(I2.2,A,I2.2,A,I2.2,A,I3.3,2(2X,I4),2X,A)')
     &   dat(5),':',dat(6),':',dat(7),'.',dat(8), kstep,kinc,label
     
      end subroutine timeStamp       
      !----------------------------------------------------------------
	  
	  
      END SUBROUTINE UEXTERNALDB
      !----------------------------------------------------------------      