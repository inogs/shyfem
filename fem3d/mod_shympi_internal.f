!
! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015    ggu     project started
!
!******************************************************************

	subroutine shympi_init_internal(my_id,n_threads)

	!use shympi

	implicit none

	integer my_id,n_threads

	include "mpif.h"

	integer ierr

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	implicit none

	include "mpif.h"

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_abort_internal

        implicit none

	include "mpif.h"

        integer ierr

	call MPI_ABORT(MPI_COMM_WORLD,ierr)

        end subroutine shympi_abort_internal

!******************************************************************

        subroutine shympi_finalize_internal

        implicit none

	include "mpif.h"

        integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)

        end subroutine shympi_finalize_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)

        implicit none

	include "mpif.h"

        integer size

	size = MPI_STATUS_SIZE

        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_internal

	implicit none

	include "mpif.h"

	integer ierr

	flush(6)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_syncronize_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_syncronize_initial

        implicit none

	include "mpif.h"

	integer my_id,nt
	integer root,ierr,i
	integer count
	integer local
	integer, allocatable :: buf(:)

        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, nt, ierr )

	allocate(buf(nt))

	root = my_id
	count = 1
	root = 0

	if( my_id == root ) then
	  do i=1,nt
	    buf(i) = i
	  end do
	end if

	call MPI_SCATTER (buf,count,MPI_INT
     +			,local,count,MPI_INT
     +			,root,MPI_COMM_WORLD,ierr)

	local = local * 2
	root = 0

	call MPI_GATHER (local,count,MPI_INT
     +			,buf,count,MPI_INT
     +			,root,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	if( my_id == root ) then
	  !write(6,*) 'mpi sync: ',nt,root,(buf(i),i=1,nt)
	  !write(6,*) 'mpi sync: ',nt,root
	end if

	deallocate(buf)

	end subroutine shympi_syncronize_initial

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_internal_i(nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi

	implicit none

	include "mpif.h"

	integer nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	integer val(nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb

        tag=1234
	ir = 0

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  call count_buffer(nlvddi,n,nc,il,g_out(:,ia),nb)
	  !write(6,*) 'ex1: ',my_id,ia,id,nc,n,nb
          call MPI_Irecv(i_buffer_out(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(3,ia)
	  call to_buffer_i(nlvddi,n,nc,il
     +		,g_in(:,ia),val,nb,i_buffer_in(:,ia))
          call MPI_Isend(i_buffer_in(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  call from_buffer_i(nlvddi,n,nc,il
     +		,g_out(:,ia),val,nb,i_buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_i

!******************************************************************

	subroutine shympi_exchange_internal_d(nlvddi,n,il,val)
	
	integer nlvddi,n
	integer il(n)
	double precision val(nlvddi,n)

        stop 'error stop shympi_exchange_internal_d: not ready'

	end subroutine shympi_exchange_internal_d

!******************************************************************

	subroutine shympi_exchange_internal_r(nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi

	implicit none

	include "mpif.h"

	integer nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	real val(nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb

        tag=1234
	ir = 0

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  call count_buffer(nlvddi,n,nc,il,g_out(:,ia),nb)
	  !write(6,*) 'ex1: ',my_id,ia,id,nc,n,nb
          call MPI_Irecv(r_buffer_out(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(3,ia)
!	  call to_buffer_r(nlvddi,n,nc,il
!     +		,g_in(:,ia),val,nb,r_buffer_in(:,ia))
          call MPI_Isend(r_buffer_in(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
!	  call from_buffer_r(nlvddi,n,nc,il
!     +		,g_out(:,ia),val,nb,r_buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_r

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_gather_i_internal(val)

	use shympi

	implicit none

        integer val

	include "mpif.h"

        integer ierr

        call MPI_GATHER (val,1,MPI_INT
     +                  ,ival,1,MPI_INT
     +                  ,0,MPI_COMM_WORLD,ierr)

        end subroutine shympi_gather_i_internal

!*******************************

        subroutine shympi_bcast_i_internal(val)

	implicit none

        integer val

	include "mpif.h"

        integer ierr

        call MPI_BCAST(val,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

        end subroutine shympi_bcast_i_internal

!******************************************************************
!******************************************************************
!******************************************************************
