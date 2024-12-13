   subroutine init3d
      use flow_data
      implicit none
      integer:: i, j, k, ierr

      call init_mesh

      call read_flow_data
      
      do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            f(i, j, k, 1) = d(i, j, k)
            f(i, j, k, 2) = d(i, j, k)*u(i, j, k)
            f(i, j, k, 3) = d(i, j, k)*v(i, j, k)
            f(i, j, k, 4) = d(i, j, k)*w(i, j, k)
          f(i, j, k, 5) = d(i, j, k)*((u(i, j, k)*u(i, j, k) + v(i, j, k)*v(i, j, k) + w(i, j, k)*w(i, j, k))*0.5d0 + Cv*T(i, j, k))
         end do
      end do
      end do
   end
