!
! implements flood fill for finite element grid
!
!*********************************************************************

	subroutine flood_fill(color)

! fills color starting from some colored points
!
! ic = color(k), k=1,nkn
!
! ic == 0 indicates a point that still has to be colored
! ic > 0  indicates a node colored with color ic
! ic < 0  indicates a point that not has to be touched
!
! at the end of the routine all points are different from 0

	use basin

	implicit none

	integer color(nkn)

	integer ngood,nbad
	integer ic
	integer it,is,isc
	integer ie,ii,k,n
	integer kn(3)
	integer coloraux(nkn)

	coloraux = 0

	nbad = count( color < 0 )

	do

	  ngood = count( color > 0 )

          write(6,*) 'influence: ',ngood,nbad,nkn,nkn-ngood-nbad

	  if( nkn-ngood-nbad == 0 ) exit

          do ie=1,nel

	    call basin_get_vertex_nodes(ie,n,kn)

            it = 0
            is = 0
            do ii=1,n
              k = kn(ii)
              if( color(k) > 0 ) then
                it = it + 1
                is = is + ii
              end if
            end do

            if( it == 0 ) then        !no node is colored
	      !nothing to do
            else if( it == 1 ) then   !one node is colored
              ic = color(kn(is))
              do ii=1,n
                coloraux(kn(ii)) = ic
              end do
            else if( it < n ) then	!two nodes are colored (in 2D)
              is = 6 - is		!this is the uncolored node
              isc = mod(is,3) + 1	!take color from this node
              k = kn(isc)
              ic = color(k)
              k = kn(is)
              coloraux(k) = ic
            end if

            if( it > 0 ) then        !check
              do ii=1,n
                k = kn(ii)
                if( color(k) == 0 .and. coloraux(k) == 0 ) then
                  write(6,*) 'internal error...... '
                  write(6,*) ie,it,is
                  stop 'error stop: internal error'
                end if
              end do
            end if

          end do

	  where( coloraux > 0 .and. color == 0 ) color = coloraux

        end do

	end

!*********************************************************************

