
  subroutine flowfield_init(flowtype)

     Use flow_data

     implicit none

     ! argument
     character(len=*), intent(in) :: flowtype

     ! local data
     integer :: i,j,k,ierr

     Istep=0

     tt=0.d0

     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
      d(i,j,k)=1.d0
      u(i,j,k)= sin(Axx(i,j,k))*cos(Ayy(i,j,k))*cos(Azz(i,j,k))
      v(i,j,k)=-cos(Axx(i,j,k))*sin(Ayy(i,j,k))*cos(Azz(i,j,k))
      w(i,j,k)=0.d0
      T(i,j,k)=1.d0
     enddo
     enddo
     enddo

     call MPI_bcast(Istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_bcast(tt, 1, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)

     if (my_id == 0) print *, "flow initialized"

  end subroutine flowfield_init