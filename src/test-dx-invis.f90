! test scheme in inviscous term
!---------------------------------------------------------------------
   subroutine NS_solver
      use flow_data
      implicit none
      real(kind=OCFD_REAL_KIND), allocatable, dimension(:) :: fpx, fmx, Scm_Hbx, hhx1, hhx2
      real(kind=OCFD_REAL_KIND):: PI, x, dx1, dx2, dx0
      integer:: i

!-----------------------------------------------------------------------
      allocate (fpx(1 - LAP:nx + LAP), fmx(1 - LAP:nx + LAP), Scm_Hbx(0:nx), hhx1(0:nx), hhx2(0:nx))

      PI = 3.14159265358979d0
      hx = 2.d0*PI/nx_global

      Scm_Hbx = 0

      do i = 1 - LAP, nx + LAP
         x = (i - 1.d0)*hx
         fpx(i) = cos(x)
         fmx(i) = cos(x)
      end do

      call OCFD2d_flux1(fpx, hhx1, nx, LAP, Scheme%Bound_index(1, 1), Scm_Hbx, Scheme%Scheme_boundary(1:2))
      call OCFD2d_flux2(fmx, hhx2, nx, LAP, Scheme%Bound_index(1, 1), Scm_Hbx, Scheme%Scheme_boundary(1:2))

      open (99, file="test-dx.dat")
      do i = 1, nx
         x = (i - 1.d0)*hx
         dx1 = (hhx1(i) - hhx1(i - 1))/hx
         dx2 = (hhx2(i) - hhx2(i - 1))/hx
         dx0 = -sin(x)
         write (99, "(I5,2E20.10)") i, dx1 - dx0, dx2 - dx0
      end do

      deallocate (fpx, fmx, Scm_Hbx, hhx1, hhx2)
   end

