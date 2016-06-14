c
c routines for non hydrostatic terms
c
c revision log :
c
c 10.05.2013    dbf     written from scratch
c 31.05.2013    dbf     written from scratch
c
c********************************************************************

	subroutine nonhydro_init

	implicit none

	integer inohyd
	real getpar

	inohyd = nint(getpar('inohyd'))
	if( inohyd /= 0 ) then
	  write(6,*) 'inohyd = ',inohyd
	  stop 'error stop nonhydro_init: cannot run non-hydrostatic'
	end if

	end

c********************************************************************

	subroutine nonhydro_get_flag(bnohyd)

	implicit none

	logical bnohyd

	bnohyd = .false.

	end

c********************************************************************

	subroutine nonhydro_adjust

	end

c********************************************************************

	subroutine nonhydro_copy

        end

c********************************************************************

	subroutine sp256wnh 

	end

c******************************************************************

	subroutine nonhydro_set_explicit 

	end 

c**********************************************************************

	subroutine nonhydro_prepare_matrix

	end

c********************************************************************

        subroutine nonhydro_adjust_value

	end

c********************************************************************

	subroutine nonhydro_correct_uveta

	end

c********************************************************************

	subroutine system_assemble_nh_3d(ie,kn,l,mat,matl,vnot)
	
	end

c********************************************************************

	subroutine system_init_nh3dmatrix

	end

c******************************************************************

        subroutine sorti(n,ra)

        end

c********************************************************************

