
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! mercury dummy routines
!
! contents :
!
! revision log :
!
! 15.05.2016	ggu	started mercury from bio3d
! 25.05.2016	ggu	changed VERS_7_5_10
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 13.06.2017	ggu	changed VERS_7_5_29
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 09.01.2020	ggu	dummy routine written
! 09.03.2020	ggu	dummy restart routines
!
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!
! notes :
!
! State variables used: (mercury) -> Donata, please adjourn
!
! State variables used: (mer
! Hg0           81      1
! Hg2           82      2
! Hg3           83      3
!
!********************************************************************

        subroutine mercury_module

! general interface to mercury module

        implicit none

        integer, save :: icall = 0

	integer imerc
	real getpar

        if( icall .le. -1 ) return

        if( icall .eq. 0 ) then
          imerc = nint(getpar('imerc'))
          if( imerc .le. 0 ) icall = -1
          if( icall .le. -1 ) return
          stop 'error stop mercury_init: mercury module is not linked'
        end if

        end

!********************************************************************

	subroutine mercury_init()
	implicit none
	end

!********************************************************************

	subroutine mercury_has_restart()
	implicit none
	end

!*************************************************************

	subroutine write_restart_mercury(iunit)
	implicit none
	integer iunit
	end

	subroutine skip_restart_mercury(iunit)
          use mercury
          use levels, only : nlvdi
          use basin, only : nkndi
        
          implicit none
          integer iunit
          integer l,k,s,flag
          integer r_ns,r_nl,r_nk,r_flag

            read(iunit) r_ns,r_nl,r_nk,r_flag
            if( r_flag/=-99991   ) goto 98
            do s=1,npstate
              read(iunit) 
            enddo
           
            read(iunit) r_ns,r_nl,r_nk,r_flag
            if(r_flag/=-99992   ) goto 98
            do s=1,nsstate
              read(iunit) 
            enddo
           
            read(iunit) r_ns,r_nl,r_nk,r_flag
            if(r_flag/=-99993   ) goto 98
            do s=1,nsolwst
              read(iunit) 
            enddo
           
            read(iunit) r_ns,r_nl,r_nk,r_flag
            if(r_flag/=-99994   ) goto 98
            do s=1,nsolsst
              read(iunit) 
            enddo
           
            read(iunit) r_ns,r_nl,r_nk,r_flag
            if(r_flag/=-99995   ) goto 98
              read(iunit) 
           
            read(iunit) r_ns,r_nl,r_nk,r_flag
            if(r_flag/=-99996   ) goto 98
              read(iunit) 
       
          return 
           
   98     continue
            write(6,*) 'error flag not conform in skip_restart_merc '
            write(6,*) r_ns,r_nl,r_nk,r_flag
            stop 'error stop skip_restart_merc'

        end

        subroutine read_restart_mercury(iunit)
          implicit none
          integer iunit
          call skip_restart_mercury(iunit)
        end

!*************************************************************

