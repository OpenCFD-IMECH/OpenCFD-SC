! Filtering, to remove high-wavenumber oscillations
! Bogey C, Bailly C,  J. Comput. Phys. 194 (2004) 194-214

     subroutine filtering(f)
        Use flow_para
        implicit none
        integer:: m, ib, ie, jb, je, kb, ke, IF_filter, Filter_X, Filter_Y, Filter_Z, Filter_scheme
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), s0, rth
        real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :, :):: f0, p
        integer, parameter:: Filter_Fo9p = 1, Filter_Fopt_shock = 2

        IF_filter = 0
        do m = 1, Para%Nfiltering
           if (mod(Istep, nint(Para%Filter(1, m))) == 0) IF_filter = 1
        end do
        if (IF_filter == 0) return    ! do not filtering in this step

        allocate (f0(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP), p(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP))
!  --------------Filtering --------------------
        if (my_id .eq. 0) print *, "filtering ......"

        do m = 1, Para%NFiltering
           if (mod(Istep, nint(Para%Filter(1, m))) == 0) then
              Filter_X = nint(Para%Filter(2, m))
              Filter_Y = nint(Para%Filter(3, m))
              Filter_Z = nint(Para%Filter(4, m))
              ib = nint(Para%Filter(5, m))
              ie = nint(Para%Filter(6, m))
              jb = nint(Para%Filter(7, m))
              je = nint(Para%Filter(8, m))
              kb = nint(Para%Filter(9, m))
              ke = nint(Para%Filter(10, m))
              ib = max(ib - i_offset(npx) + 1, 1)               ! transform global index to local index
              ie = min(ie - i_offset(npx) + 1, nx)
              jb = max(jb - j_offset(npy) + 1, 1)
              je = min(je - j_offset(npy) + 1, ny)
              kb = max(kb - k_offset(npz) + 1, 1)
              ke = min(ke - k_offset(npz) + 1, nz)

              Filter_scheme = nint(Para%Filter(11, m))     ! 1 9-point filtering Fo9P;   2 shock-capturing filtering (more robust but more dissipative)
              s0 = Para%Filter(12, m)               ! filtering amplitude (s=1 for default)
              rth = Para%Filter(13, m)              ! flitering threshold (1E-5 for default);   smaller: more dissipation (robust) rth=0: 100% filtering

              if (Filter_X == 1) then
                 if (Filter_scheme == Filter_Fo9p) then
                    call filter_x3d(f, f0, s0, ib, ie, jb, je, kb, ke)
                 else if (Filter_scheme == Filter_Fopt_shock) then
                    call filter_x3d_shock(f, f0, p, s0, rth, ib, ie, jb, je, kb, ke)
                 end if
              end if

              if (Filter_Y == 1) then
                 if (Filter_scheme == Filter_Fo9p) then
                    call filter_y3d(f, f0, s0, ib, ie, jb, je, kb, ke)
                 else if (Filter_scheme == Filter_Fopt_shock) then
                    call filter_y3d_shock(f, f0, p, s0, rth, ib, ie, jb, je, kb, ke)
                 end if
              end if

              if (Filter_Z == 1) then
                 if (Filter_scheme == Filter_Fo9p) then
                    call filter_z3d(f, f0, s0, ib, ie, jb, je, kb, ke)
                 else if (Filter_scheme == Filter_Fopt_shock) then
                    call filter_z3d_shock(f, f0, p, s0, rth, ib, ie, jb, je, kb, ke)
                 end if
              end if
           end if
        end do
        deallocate (f0, p)
     end

!---------------------------------------------------
     subroutine filter_x3d(f, f0, s0, ib, ie, jb, je, kb, ke)
        Use flow_para
        implicit none

        integer:: i, j, k, m, ib, ie, jb, je, kb, ke, ib1, ie1
        real(kind=OCFD_REAL_KIND):: f0(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), s0
        real(kind=OCFD_REAL_KIND), parameter::  d0 = 0.243527493120d0, d1 = -0.204788880640d0, d2 = 0.120007591680d0, &
                                               d3 = -0.045211119360d0, d4 = 0.008228661760d0

        ib1 = ib; ie1 = ie
        if (npx == 0 .and. Para%Iperiodic_X .ne. 1) ib1 = max(ib, 6)
        if (npx == npx0 - 1 .and. Para%Iperiodic_X .ne. 1) ie1 = min(ie, nx - 5)

        do m = 1, Nvars
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f0(i, j, k) = f(i, j, k, m)
              end do
              end do
           end do

           call exchange_boundary_x(f0)

           do k = kb, ke
              do j = jb, je
              do i = ib1, ie1
f(i, j, k, m) = f0(i, j, k) - s0*(d0*f0(i, j, k) + d1*(f0(i - 1, j, k) + f0(i + 1, j, k)) + d2*(f0(i - 2, j, k) + f0(i + 2, j, k)) &
                                                  + d3*(f0(i - 3, j, k) + f0(i + 3, j, k)) + d4*(f0(i - 4, j, k) + f0(i + 4, j, k)))
              end do
              end do
           end do
        end do

     end

!---------------------------------------------------
     subroutine filter_y3d(f, f0, s0, ib, ie, jb, je, kb, ke)
        Use flow_para
        implicit none

        integer:: i, j, k, m, ib, ie, jb, je, kb, ke, jb1, je1
        real(kind=OCFD_REAL_KIND):: f0(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), s0
        real(kind=OCFD_REAL_KIND), parameter::  d0 = 0.243527493120d0, d1 = -0.204788880640d0, d2 = 0.120007591680d0, &
                                               d3 = -0.045211119360d0, d4 = 0.008228661760d0

        jb1 = jb; je1 = je
        if (npy == 0 .and. Para%Iperiodic_Y .ne. 1) jb1 = max(jb, 6)
        if (npy == npy0 - 1 .and. Para%Iperiodic_Y .ne. 1) je1 = min(je, ny - 5)

        do m = 1, Nvars
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f0(i, j, k) = f(i, j, k, m)
              end do
              end do
           end do

           call exchange_boundary_y(f0)

           do k = kb, ke
              do j = jb1, je1
              do i = ib, ie
f(i, j, k, m) = f0(i, j, k) - s0*(d0*f0(i, j, k) + d1*(f0(i, j - 1, k) + f0(i, j + 1, k)) + d2*(f0(i, j - 2, k) + f0(i, j + 2, k)) &
                                                  + d3*(f0(i, j - 3, k) + f0(i, j + 3, k)) + d4*(f0(i, j - 4, k) + f0(i, j + 4, k)))
              end do
              end do
           end do
        end do

     end

!---------------------------------------------------
     subroutine filter_z3d(f, f0, s0, ib, ie, jb, je, kb, ke)
        Use flow_para
        implicit none

        integer:: i, j, k, m, ib, ie, jb, je, kb, ke, kb1, ke1
        real(kind=OCFD_REAL_KIND):: f0(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), s0
        real(kind=OCFD_REAL_KIND), parameter::  d0 = 0.243527493120d0, d1 = -0.204788880640d0, d2 = 0.120007591680d0, &
                                               d3 = -0.045211119360d0, d4 = 0.008228661760d0
        kb1 = kb; ke1 = ke
        if (npz == 0 .and. Para%Iperiodic_Z .ne. 1) kb1 = max(kb, 6)
        if (npz == npz0 - 1 .and. Para%Iperiodic_Z .ne. 1) ke1 = min(ke, nz - 5)

        do m = 1, Nvars
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f0(i, j, k) = f(i, j, k, m)
              end do
              end do
           end do

           call exchange_boundary_z(f0)

           do k = kb1, ke1
              do j = jb, je
              do i = ib, ie
f(i, j, k, m) = f0(i, j, k) - s0*(d0*f0(i, j, k) + d1*(f0(i, j, k - 1) + f0(i, j, k + 1)) + d2*(f0(i, j, k - 2) + f0(i, j, k + 2)) &
                                                  + d3*(f0(i, j, k - 3) + f0(i, j, k + 3)) + d4*(f0(i, j, k - 4) + f0(i, j, k + 4)))
              end do
              end do
           end do
        end do

     end

!------------------------------------------------------------
! Shock cpaturing filtering

     subroutine filter_x3d_shock(f, f0, p, s0, rth, ib, ie, jb, je, kb, ke)
        Use flow_para
        implicit none

        integer:: i, j, k, m, ib, ie, jb, je, kb, ke, ib1, ie1
        real(kind=OCFD_REAL_KIND), dimension(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP):: f0, p
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), sc0(0:nx + 1), s0, rth, dp0, dp1, dp2, ri, sc1, sc2
        real(kind=OCFD_REAL_KIND), parameter::   c1 = -0.210383d0, c2 = 0.039617d0

        ib1 = ib; ie1 = ie
        if (npx == 0 .and. Para%Iperiodic_X .ne. 1) ib1 = max(ib, 4)
        if (npx == npx0 - 1 .and. Para%Iperiodic_X .ne. 1) ie1 = min(ie, nx - 3)

        do k = 1, nz
           do j = 1, ny
           do i = 1, nx
     p(i, j, k) = (f(i, j, k, 5) - 0.5d0*(f(i, j, k, 2)**2 + f(i, j, k, 3)**2 + f(i, j, k, 4)**2)/f(i, j, k, 1))*(Para%gamma - 1.d0)
           end do
           end do
        end do

        call exchange_boundary_x(p)

        do m = 1, Nvars
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f0(i, j, k) = f(i, j, k, m)
              end do
              end do
           end do

           call exchange_boundary_x(f0)

           do k = kb, ke
              do j = jb, je

                 do i = ib1 - 1, ie1 + 1
                    dp0 = 0.25d0*(-p(i + 1, j, k) + 2.d0*p(i, j, k) - p(i - 1, j, k))
                    dp1 = 0.25d0*(-p(i + 2, j, k) + 2.d0*p(i + 1, j, k) - p(i, j, k))
                    dp2 = 0.25d0*(-p(i, j, k) + 2.d0*p(i - 1, j, k) - p(i - 2, j, k))
                    ri = 0.5d0*((dp0 - dp1)**2 + (dp0 - dp2)**2)/(p(i, j, k)*p(i, j, k)) + 1.d-16
                    sc0(i) = 0.5d0*(1.d0 - rth/ri + abs(1.d0 - rth/ri))
                 end do

                 do i = ib1, ie1
                    sc1 = 0.5d0*(sc0(i) + sc0(i + 1))
                    sc2 = 0.5d0*(sc0(i) + sc0(i - 1))

               f(i, j, k, m) = f0(i, j, k) - s0*(Sc1*(c1*(f0(i + 1, j, k) - f0(i, j, k)) + c2*(f0(i + 2, j, k) - f0(i - 1, j, k))) &
                                                - Sc2*(c1*(f0(i, j, k) - f0(i - 1, j, k)) + c2*(f0(i + 1, j, k) - f0(i - 2, j, k))))
                 end do
              end do
           end do
        end do

     end

!---------------------------------------------------
     subroutine filter_y3d_shock(f, f0, p, s0, rth, ib, ie, jb, je, kb, ke)
        Use flow_para
        implicit none

        integer:: i, j, k, m, ib, ie, jb, je, kb, ke, jb1, je1
        real(kind=OCFD_REAL_KIND), dimension(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP):: f0, p
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), sc0(0:ny + 1), s0, rth, dp0, dp1, dp2, ri, sc1, sc2
        real(kind=OCFD_REAL_KIND), parameter::   c1 = -0.210383d0, c2 = 0.039617d0

        jb1 = jb; je1 = je
        if (npy == 0 .and. Para%Iperiodic_Y .ne. 1) jb1 = max(jb, 4)
        if (npy == npy0 - 1 .and. Para%Iperiodic_Y .ne. 1) je1 = min(je, ny - 3)

        do k = 1, nz
           do j = 1, ny
           do i = 1, nx
     p(i, j, k) = (f(i, j, k, 5) - 0.5d0*(f(i, j, k, 2)**2 + f(i, j, k, 3)**2 + f(i, j, k, 4)**2)/f(i, j, k, 1))*(Para%gamma - 1.d0)
           end do
           end do
        end do

        call exchange_boundary_y(p)

        do m = 1, Nvars
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f0(i, j, k) = f(i, j, k, m)
              end do
              end do
           end do

           call exchange_boundary_y(f0)

           do k = kb, ke
              do i = ib, ie

                 do j = jb1 - 1, je1 + 1
                    dp0 = 0.25d0*(-p(i, j + 1, k) + 2.d0*p(i, j, k) - p(i, j - 1, k))
                    dp1 = 0.25d0*(-p(i, j + 2, k) + 2.d0*p(i, j + 1, k) - p(i, j, k))
                    dp2 = 0.25d0*(-p(i, j, k) + 2.d0*p(i, j - 1, k) - p(i, j - 2, k))
                    ri = 0.5d0*((dp0 - dp1)**2 + (dp0 - dp2)**2)/(p(i, j, k)*p(i, j, k)) + 1.d-16
                    sc0(j) = 0.5d0*(1.d0 - rth/ri + abs(1.d0 - rth/ri))
                 end do

                 do j = jb1, je1
                    sc1 = 0.5d0*(sc0(j) + sc0(j + 1))
                    sc2 = 0.5d0*(sc0(j) + sc0(j - 1))

               f(i, j, k, m) = f0(i, j, k) - s0*(Sc1*(c1*(f0(i, j + 1, k) - f0(i, j, k)) + c2*(f0(i, j + 2, k) - f0(i, j - 1, k))) &
                                                - Sc2*(c1*(f0(i, j, k) - f0(i, j - 1, k)) + c2*(f0(i, j + 1, k) - f0(i, j - 2, k))))
                 end do

              end do
           end do
        end do

     end
!---------------------------------------------------
     subroutine filter_z3d_shock(f, f0, p, s0, rth, ib, ie, jb, je, kb, ke)
        Use flow_para
        implicit none

        integer:: i, j, k, m, ib, ie, jb, je, kb, ke, kb1, ke1
        real(kind=OCFD_REAL_KIND), dimension(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP):: f0, p
        real(kind=OCFD_REAL_KIND):: f(nx, ny, nz, Nvars), sc0(0:nz + 1), s0, rth, dp0, dp1, dp2, ri, sc1, sc2
        real(kind=OCFD_REAL_KIND), parameter::   c1 = -0.210383d0, c2 = 0.039617d0

        kb1 = kb; ke1 = ke
        if (npz == 0 .and. Para%Iperiodic_Z .ne. 1) kb1 = max(kb, 4)
        if (npz == npz0 - 1 .and. Para%Iperiodic_Z .ne. 1) ke1 = min(ke, nz - 3)

        do k = 1, nz
           do j = 1, ny
           do i = 1, nx
     p(i, j, k) = (f(i, j, k, 5) - 0.5d0*(f(i, j, k, 2)**2 + f(i, j, k, 3)**2 + f(i, j, k, 4)**2)/f(i, j, k, 1))*(Para%gamma - 1.d0)
           end do
           end do
        end do

        call exchange_boundary_z(p)

        do m = 1, Nvars
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f0(i, j, k) = f(i, j, k, m)
              end do
              end do
           end do

           call exchange_boundary_z(f0)

           do j = jb, je
           do i = ib, ie

              do k = kb1 - 1, ke1 + 1
                 dp0 = 0.25d0*(-p(i, j, k + 1) + 2.d0*p(i, j, k) - p(i, j, k - 1))
                 dp1 = 0.25d0*(-p(i, j, k + 2) + 2.d0*p(i, j, k + 1) - p(i, j, k))
                 dp2 = 0.25d0*(-p(i, j, k) + 2.d0*p(i, j, k - 1) - p(i, j, k - 2))
                 ri = 0.5d0*((dp0 - dp1)**2 + (dp0 - dp2)**2)/(p(i, j, k)*p(i, j, k)) + 1.d-16
                 sc0(k) = 0.5d0*(1.d0 - rth/ri + abs(1.d0 - rth/ri))
              end do

              do k = kb1, ke1
                 sc1 = 0.5d0*(sc0(k) + sc0(k + 1))
                 sc2 = 0.5d0*(sc0(k) + sc0(k - 1))

               f(i, j, k, m) = f0(i, j, k) - s0*(Sc1*(c1*(f0(i, j, k + 1) - f0(i, j, k)) + c2*(f0(i, j, k + 2) - f0(i, j, k - 1))) &
                                                - Sc2*(c1*(f0(i, j, k) - f0(i, j, k - 1)) + c2*(f0(i, j, k + 1) - f0(i, j, k - 2))))
              end do

           end do
           end do
        end do

     end

!------------------------------------------------------------
