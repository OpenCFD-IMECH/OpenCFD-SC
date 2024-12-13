! 3D Compressible N-S Finite Difference Solver
! CopyRight by Li Xinliang  Email: lixl@imech.ac.cn
! Version 2.x   2021-2
!---------------------------------------------------------------------
   subroutine NS_solver
      use flow_data
      implicit none
      real*8:: wall_time,wall_time_beg,wall_time_end
      integer:: KRK, i, j, k, m, ierr
!-----------------------------------------------------------------------
      call allocate_flow_data        ! f, fn, d,u,v,w,T, Axx,Ayy,....,Ajac
      call allocate_inviscous_data   ! work data for inviscous terms
      call allocate_vicous_data      ! work data for viscous terms
      call init3d
      call OCFD_bc          ! boundary condition

      if (my_id .eq. 0) print *, 'init ok'
      wall_time = MPI_wtime()
      wall_time_beg=wall_time
!c-----------------------------------------------------------------------
      do while (tt < Para%End_time - Para%dt*1.d-4)       ! -dt*1.d-4 , considering rounding error

         do m = 1, 5
            do k = 1, nz
               do j = 1, ny
               do i = 1, nx
                  fn(i, j, k, m) = f(i, j, k, m)
               end do
               end do
            end do
         end do

         do KRK = 1, 3            ! 3-step Runge-Kutta

            call exchange_boundary_xyz(d)
            call exchange_boundary_xyz(u)
            call exchange_boundary_xyz(v)
            call exchange_boundary_xyz(w)
            call exchange_boundary_xyz(T)

            if (KRK .eq. 1 .and. Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
               call comput_Rhybrid           ! comput only KRK=1
            end if

            if (Para%IF_Scheme_Character .eq. 0) then
               call du_inviscous    ! inviscous term  (non-character type)
            else
               call du_inviscous_Character     ! inviscous term (character type)
            end if

            if (Para%IF_Viscous .eq. 1) then
               call du_viscous
            end if

            call OCFD_time_adv_RK3(KRK)          ! time advance (3-step RK)

            call comput_duvwT           ! comput d,u,v,w,T
            call OCFD_bc             ! boundary condition

         end do

!c ---------------4. Loop of t ------------------

         Istep = Istep + 1
         tt = tt + Para%dt

         call filtering(f)

         if (mod(Istep, Para%Istep_show) .eq. 0) then
            call show_flow_msg(wall_time)  ! CPU time, Total energy, Kinetic energy, Total entropy
         end if

         call OCFD_analysis
         if (mod(Istep, Para%Istep_Save) .eq. 0) call save_flow_data

!       call MPI_barrier(MPI_COMM_WORLD,ierr)

      end do
      wall_time_end = MPI_wtime()
      if (my_id .eq. 0)  print*,'total wtime cost:',wall_time_end-wall_time_beg
!c---------------------------------------------------------------------
      call save_flow_data
      if (my_id .eq. 0) print *, 'OK The END of opencfd'

      call deallocate_flow_data
      call deallocate_inviscous_data
      call deallocate_vicous_data
   end

!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
   subroutine allocate_flow_data
      use flow_data
      implicit none
!-----------------------------------------------------------------------

! allocate flow data space
      allocate (f(nx, ny, nz, 5), fn(nx, ny, nz, 5), &
                du(nx, ny, nz, 5), Amu(nx, ny, nz))

      allocate (d(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                u(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                v(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                w(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                T(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP))

      d = 1.d0; u = 1.d0; v = 1.d0; w = 1.d0; T = 1.d0       ! initial as 1.0
      Amu = 0.d0                                     ! initial as 0
      du = 0.d0                       ! initial as 0
!        Amu_t=0.d0

!-----Coordinate and Jacobian coefficients   (Akx1=Akx/Ajac) -----------
      allocate (Axx(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Ayy(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Azz(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Akx(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aky(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Akz(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aix(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aiy(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aiz(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Asx(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Asy(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Asz(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Ajac(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Akx1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aky1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Akz1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aix1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aiy1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Aiz1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Asx1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Asy1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), &
                Asz1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP))

      if (Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
         allocate (Rhybrid(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP))
         Rhybrid = 0.d0
      end if

!  ------initial as 1.0 --------
      Axx = 1.d0; Ayy = 1.d0; Azz = 1.d0
      Akx = 1.d0; Aky = 1.d0; Akz = 1.d0
      Aix = 1.d0; Aiy = 1.d0; Aiz = 1.d0
      Asx = 1.d0; Asy = 1.d0; Asz = 1.d0
      Ajac = 1.d0
      Akx1 = 1.d0; Aky1 = 1.d0; Akz1 = 1.d0
      Aix1 = 1.d0; Aiy1 = 1.d0; Aiz1 = 1.d0
      Asx1 = 1.d0; Asy1 = 1.d0; Asz1 = 1.d0
   end

   subroutine deallocate_flow_data
      use flow_data
      implicit none
      deallocate (f, fn, du, Amu, d, u, v, T, Axx, Ayy, Akx, Aky, Aix, Aiy, Ajac, Akx1, Aky1, Aix1, Aiy1)
      if (Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
         deallocate (Rhybrid)
      end if
   end
