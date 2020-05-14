
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014,2019  Georg Umgiesser
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

	!integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	!common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

! revision log :
!
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 16.02.2019	ggu	changed VERS_7_5_60

	integer nkn,nel,ngr,mbw
	common /nbasin/ nkn,nel,ngr,mbw
	save /nbasin/

	!integer nbc

