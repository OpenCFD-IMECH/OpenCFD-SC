! Post-analysis codes
!-----------------------------------------------------------------
 subroutine OCFD_analysis
!-------------------------------------
    Use flow_data
    implicit none
    integer m, Kana, Kstep_ana
    real(kind=OCFD_REAL_KIND):: ana_data(50)
    do m = 1, Para%ANA_Number
       Kana = nint(Para%ANA_Para(1, m))
       Kstep_ana = nint(Para%ANA_Para(2, m))
       ana_data(:) = Para%ANA_Para(3:52, m)            ! 50 elements in ana_data

       if (mod(Istep, Kstep_ana) .eq. 0) then
          select case (Kana)
          case (OCFD_ANA_time_average)
             call ana_time_average(ana_data)
          case (OCFD_ANA_Q)
             call ana_getQ(ana_data)
          case (OCFD_ANA_BOX)
             call ana_box
          case (OCFD_ANA_SAVEDATA)
             call ana_savedata(ana_data)
          case (OCFD_ANA_Corner)
             call Ana_corner
          case (OCFD_ANA_USER)
             call ana_user(ana_data)             ! user defined analysis code
          case default
             if (my_id .eq. 0) print *, "This analysis code is not supported!"
          end select
       end if
    end do
 end
!---------------------------------------------------------------

!-----------------------------------------------------------------
! code for Time averaging-----------------------------------------------
 subroutine ana_time_average(ana_data)
    Use flow_data
    implicit none
    integer i, j, k, ierr, Kstep_save_average
    real(kind=OCFD_REAL_KIND):: ana_data(50)
    real(kind=OCFD_REAL_KIND), save, allocatable, dimension(:, :, :):: dm, um, vm, wm, Tm
    real(kind=OCFD_REAL_KIND), save:: t_average = 0
    integer, save:: K_average, Iflag_average = 0
    logical ex
    character(len=50) filename
!----------------------------------------------------------------
    Kstep_save_average = nint(ana_data(1))

    if (Iflag_average .eq. 0) then   ! run only in the first time
       Iflag_average = 1
       allocate (dm(nx, ny, nz), um(nx, ny, nz), vm(nx, ny, nz), wm(nx, ny, nz), Tm(nx, ny, nz))
       if (my_id .eq. 0) inquire (file="opencfd.dat.average", exist=ex)
       call MPI_bcast(ex, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
       if (ex) then

          if (my_id .eq. 0) then
             print *, "opencfd.dat.average is found,  read it ......"
             open (88, file="opencfd.dat.average", form="unformatted")
             read (88) K_average, t_average
             print *, "K_average=", K_average, "t_average=", t_average
          end if
          call MPI_bcast(K_average, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          call MPI_bcast(t_average, 1, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)
          call read_3d(88, dm, nx, ny, nz, nx_global, ny_global, nz_global)
          call read_3d(88, um, nx, ny, nz, nx_global, ny_global, nz_global)
          call read_3d(88, vm, nx, ny, nz, nx_global, ny_global, nz_global)
          call read_3d(88, wm, nx, ny, nz, nx_global, ny_global, nz_global)
          call read_3d(88, Tm, nx, ny, nz, nx_global, ny_global, nz_global)
          if (my_id .eq. 0) close (88)

       else
          if (my_id .eq. 0) print *, "can not find opencfd.dat.average, initial the average as zero......"
          t_average = 0.d0
          K_average = 0
          dm = 0.d0
          um = 0.d0
          vm = 0.d0
          wm = 0.d0
          Tm = 0.d0
       end if
    end if
!-------------------------------------------------------
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       dm(i, j, k) = (K_average*dm(i, j, k) + d(i, j, k))/(K_average + 1.d0)
       um(i, j, k) = (K_average*um(i, j, k) + u(i, j, k))/(K_average + 1.d0)
       vm(i, j, k) = (K_average*vm(i, j, k) + v(i, j, k))/(K_average + 1.d0)
       wm(i, j, k) = (K_average*wm(i, j, k) + w(i, j, k))/(K_average + 1.d0)
       Tm(i, j, k) = (K_average*Tm(i, j, k) + T(i, j, k))/(K_average + 1.d0)
    end do
    end do
    end do

    K_average = K_average + 1
    t_average = t_average + Para%dt
    if (my_id .eq. 0) print *, "average flow ...", k_average

!---------------------------------------------------------
    if (mod(Istep, Kstep_save_average) .eq. 0) then
       if (my_id .eq. 0) then
          write (filename, "('OCFD'I8.8'.dat.average')") Istep
          print *, "write average file ...", filename
          open (99, file=filename, form='unformatted')
          write (99) K_average, t_average
       end if

       call write_3d(99, dm, nx, ny, nz, nx_global, ny_global, nz_global)
       call write_3d(99, um, nx, ny, nz, nx_global, ny_global, nz_global)
       call write_3d(99, vm, nx, ny, nz, nx_global, ny_global, nz_global)
       call write_3d(99, wm, nx, ny, nz, nx_global, ny_global, nz_global)
       call write_3d(99, Tm, nx, ny, nz, nx_global, ny_global, nz_global)

       if (my_id .eq. 0) then
          close (99)
          print *, 'write average data ok'
       end if
    end if
 end

! Compute Q (the second invarient of velocity grident) and Lamda2 (the medimum eigenvalue of Sij*Sij+Wij*Wij)
 subroutine ana_getQ(ana_data)
    use flow_data
    use viscous_data                ! work data uk,vk,wk,ui,vi,wi,us,vs,ws,TK
    implicit none
    integer i, j, k
    real(kind=OCFD_REAL_KIND)::  ana_data(50)
    real(kind=OCFD_REAL_KIND)::  ux, uy, uz, vx, vy, vz, wx, wy, wz
    character(len=100) filename
!------------------------------------------------------------------------
    call exchange_boundary_xyz(u)
    call exchange_boundary_xyz(v)
    call exchange_boundary_xyz(w)

    call OCFD_dx0(u, uk, Scheme%Scheme_Vis)
    call OCFD_dx0(v, vk, Scheme%Scheme_Vis)
    call OCFD_dx0(w, wk, Scheme%Scheme_Vis)
    call OCFD_dy0(u, ui, Scheme%Scheme_Vis)
    call OCFD_dy0(v, vi, Scheme%Scheme_Vis)
    call OCFD_dy0(w, wi, Scheme%Scheme_Vis)
    call OCFD_dz0(u, us, Scheme%Scheme_Vis)
    call OCFD_dz0(v, vs, Scheme%Scheme_Vis)
    call OCFD_dz0(w, ws, Scheme%Scheme_Vis)
!c-----------------------------------
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       ux = uk(i, j, k)*Akx(i, j, k) + ui(i, j, k)*Aix(i, j, k) + us(i, j, k)*Asx(i, j, k)
       vx = vk(i, j, k)*Akx(i, j, k) + vi(i, j, k)*Aix(i, j, k) + vs(i, j, k)*Asx(i, j, k)
       wx = wk(i, j, k)*Akx(i, j, k) + wi(i, j, k)*Aix(i, j, k) + ws(i, j, k)*Asx(i, j, k)
       uy = uk(i, j, k)*Aky(i, j, k) + ui(i, j, k)*Aiy(i, j, k) + us(i, j, k)*Asy(i, j, k)
       vy = vk(i, j, k)*Aky(i, j, k) + vi(i, j, k)*Aiy(i, j, k) + vs(i, j, k)*Asy(i, j, k)
       wy = wk(i, j, k)*Aky(i, j, k) + wi(i, j, k)*Aiy(i, j, k) + ws(i, j, k)*Asy(i, j, k)
       uz = uk(i, j, k)*Akz(i, j, k) + ui(i, j, k)*Aiz(i, j, k) + us(i, j, k)*Asz(i, j, k)
       vz = vk(i, j, k)*Akz(i, j, k) + vi(i, j, k)*Aiz(i, j, k) + vs(i, j, k)*Asz(i, j, k)
       wz = wk(i, j, k)*Akz(i, j, k) + wi(i, j, k)*Aiz(i, j, k) + ws(i, j, k)*Asz(i, j, k)
       TK(i, j, k) = ux*vy + ux*wz + vy*wz - uy*vx - uz*wx - vz*wy  !! TK=Q=II(UX)
    end do
    end do
    end do
!-----save data---------------------------
    if (my_id .eq. 0) then
       write (filename, "('Q-'I7.7'.dat')") Istep
       open (106, file=filename, form="unformatted")
    end if
    call write_3d(106, TK, nx, ny, nz, nx_global, ny_global, nz_global)
    if (my_id .eq. 0) close (106)
!-----------------------------------
 end

!--------------------------------------------------------------------------------------------
!------------------------------------------------
! Copyright by Li Xinliang
! Compute the statistics data of isotropic turbulence, such as Re_lamda, Mt, Skewness and Flatness factors of u, ux ......
! Ref: JCP 13(5):1415,2001
 subroutine ana_box
    use flow_data
    use viscous_data                ! work data uk,vk,wk,ui,vi,wi,us,vs,ws,TK
    implicit none

    integer i, j, k, ierr
    real(kind=OCFD_REAL_KIND):: ux,uy,uz,vx,vy,vz,wx,wy,wz,amu_aver, c_aver, d_aver, T_aver, u_rms, v_rms, w_rms,ens,&
                                VV_rms, Ux_rms, d_rms, E_kinetic, ux_skewness, ux_flatness, u_skewness, u_flatness, &
                                Alamda, Re_lamda, Alamdax, Re_lamdax, Amt, atmp1(100), atmp0(100)

!---------------------------------------------------------
! computing statistical data ....
    c_aver = 0.d0; d_aver = 0.d0; Amu_aver = 0.d0; T_aver = 0.d0    ! Mean sound speed, density, viscous, Temperature
    u_rms = 0.d0; v_rms = 0.d0; w_rms = 0.d0; ux_rms = 0.d0      ! RMS of velocity fluctuations
    ux_skewness = 0.d0; ux_flatness = 0.d0; u_skewness = 0.d0; u_flatness = 0.d0  ! Skewness & Flatness factor of velocities
    E_kinetic = 0.d0

    ens=0.d0

    call exchange_boundary_xyz(u)
    call OCFD_dx0(u, uk, Scheme%Scheme_Vis)
    call OCFD_dy0(u, ui, Scheme%Scheme_Vis)
    call OCFD_dz0(u, us, Scheme%Scheme_Vis)

    call exchange_boundary_xyz(v)
    call OCFD_dx0(v, vk, Scheme%Scheme_Vis)
    call OCFD_dy0(v, vi, Scheme%Scheme_Vis)
    call OCFD_dz0(v, vs, Scheme%Scheme_Vis)

    call exchange_boundary_xyz(w)
    call OCFD_dx0(w, wk, Scheme%Scheme_Vis)
    call OCFD_dy0(w, wi, Scheme%Scheme_Vis)
    call OCFD_dz0(w, ws, Scheme%Scheme_Vis)

    call comput_Amu

! Statistics
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       ux = uk(i, j, k)*Akx(i, j, k) + ui(i, j, k)*Aix(i, j, k) + us(i, j, k)*Asx(i, j, k)
       uy = uk(i, j, k)*Aky(i, j, k) + ui(i, j, k)*Aiy(i, j, k) + us(i, j, k)*Asy(i, j, k)
       uz = uk(i, j, k)*Akz(i, j, k) + ui(i, j, k)*Aiz(i, j, k) + us(i, j, k)*Asz(i, j, k)
       vx = vk(i, j, k)*Akx(i, j, k) + vi(i, j, k)*Aix(i, j, k) + vs(i, j, k)*Asx(i, j, k)
       vy = vk(i, j, k)*Aky(i, j, k) + vi(i, j, k)*Aiy(i, j, k) + vs(i, j, k)*Asy(i, j, k)
       vz = vk(i, j, k)*Akz(i, j, k) + vi(i, j, k)*Aiz(i, j, k) + vs(i, j, k)*Asz(i, j, k)
       wx = wk(i, j, k)*Akx(i, j, k) + wi(i, j, k)*Aix(i, j, k) + ws(i, j, k)*Asx(i, j, k)
       wy = wk(i, j, k)*Aky(i, j, k) + wi(i, j, k)*Aiy(i, j, k) + ws(i, j, k)*Asy(i, j, k)
       wz = wk(i, j, k)*Akz(i, j, k) + wi(i, j, k)*Aiz(i, j, k) + ws(i, j, k)*Asz(i, j, k)
       c_aver = c_aver + sqrt(T(i, j, k))/Para%Ma
       d_aver = d_aver + d(i, j, k)
       Amu_aver = Amu_aver + Amu(i, j, k)
       T_aver = T_aver + T(i, j, k)
       u_rms = u_rms + u(i, j, k)*u(i, j, k)
       v_rms = v_rms + v(i, j, k)*v(i, j, k)
       w_rms = w_rms + w(i, j, k)*w(i, j, k)

       Ux_rms = Ux_rms + ux*ux
       Ux_skewness = Ux_skewness + ux**3
       Ux_flatness = Ux_flatness + ux**4
       u_skewness = u_skewness + u(i, j, k)**3
       u_flatness = u_flatness + u(i, j, k)**4
       E_kinetic = E_kinetic + (u(i, j, k)*u(i, j, k) + v(i, j, k)*v(i, j, k) + w(i, j, k)*w(i, j, k))*d(i, j, k)*0.5d0

       ens = ens + (wy-vz)**2+(uz-wx)**2+(vx-uy)**2
    end do
    end do
    end do

    atmp1(1) = u_rms; atmp1(2) = v_rms; atmp1(3) = w_rms; atmp1(4) = Ux_rms
    atmp1(5) = Ux_skewness; atmp1(6) = Ux_flatness
    atmp1(7) = d_aver; atmp1(8) = c_aver; atmp1(9) = Amu_aver
    atmp1(10) = E_kinetic; atmp1(11) = T_aver; atmp1(12) = u_skewness; atmp1(13) = u_flatness
    atmp1(14) = ens

    call MPI_ALLREDUCE(atmp1, atmp0, 14, OCFD_DATA_TYPE, MPI_SUM, MPI_Comm_World, ierr)

    atmp0 = atmp0/(nx_global*ny_global*nz_global)

    u_rms = atmp0(1); v_rms = atmp0(2); w_rms = atmp0(3); Ux_rms = atmp0(4)
    Ux_skewness = atmp0(5); Ux_flatness = atmp0(6); d_aver = atmp0(7)
    c_aver = atmp0(8); Amu_aver = atmp0(9); E_kinetic = atmp0(10)
    T_aver = atmp0(11); u_skewness = atmp0(12); u_flatness = atmp0(13)
    ens = 0.5d0*atmp0(14)

!-----------------------------------------
    VV_rms = sqrt((u_rms + v_rms + w_rms)/3.d0)  !  RMS of velocity fluctures
    Amt = sqrt(u_rms + v_rms + w_rms)/c_aver     ! Turbulent Mach number (see JCP 13(5):1416)  (need not div by 3)

    u_rms = sqrt(u_rms); v_rms = sqrt(v_rms); w_rms = sqrt(w_rms); Ux_rms = sqrt(Ux_rms)

    Alamda = VV_rms/Ux_rms   ! Taylor scale, the same as Samtaney et al.
    Alamdax = u_rms/Ux_rms

    Re_lamda = VV_rms*Alamda*d_aver/Amu_aver
    Re_lamdax = VV_rms*Alamdax*d_aver/Amu_aver

    Ux_skewness = Ux_skewness/Ux_rms**3
    Ux_flatness = Ux_flatness/Ux_rms**4
    U_skewness = U_skewness/U_rms**3
    U_flatness = u_flatness/U_rms**4
!--------d_rms----------------------------------------------

    d_rms = 0.d0
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       d_rms = d_rms + (d(i, j, k) - d_aver)**2
    end do
    end do
    end do
    call MPI_ALLREDUCE(d_rms, atmp1(1), 1, OCFD_DATA_TYPE, MPI_SUM, MPI_Comm_World, ierr)
    d_rms = sqrt(atmp1(1)/(nx_global*ny_global*nz_global))
!---------------------------------------------------------
    if (my_id .eq. 0) then

       write(*,'(4(1X,A13))')"tt", "Re_lamda", "Re_lamdax", "Amt"
       write(*,'(4(1X,E13.6E2))')tt, Re_lamda, Re_lamdax, Amt
       write(*,'(4(1X,A13))')"vv_rms", "u_rms", "v_rms", "w_rms"
       write(*,'(4(1X,E13.6E2))')vv_rms, u_rms, v_rms, w_rms
       write(*,'(4(1X,A13))')"u_skewness", "u_flatness", "ux_skewness", "ux_flatness"
       write(*,'(4(1X,E13.6E2))')u_skewness, u_flatness, ux_skewness, ux_flatness
       write(*,'(4(1X,A13))')"E_kinetic","Enstrophy", "d_rms", "T_aver"
       write(*,'(4(1X,E13.6E2))')E_kinetic,ens,d_rms,T_aver
       write(*,'(4(1X,A13))')"Amu_aver", "d_aver", "Ux_rms", "Alamdax"
       write(*,'(4(1X,E13.6E2))')Amu_aver, d_aver, Ux_rms, Alamdax

       open (33, file='statistics.dat', position='append')
       write (33, '(20(E20.13E2))') tt, Re_lamda, Re_lamdax, Amt, vv_rms, u_rms, v_rms, w_rms, &
          u_skewness, u_flatness, ux_skewness, ux_flatness, &
          E_kinetic,ens,d_rms, T_aver, Amu_aver, d_aver, Ux_rms, Alamdax
       close (33)
    end if
 end

!-------------------------------------------------------------------------------
! Post analysis code for Compression corner
! save pw;  u, T at j=2
 subroutine Ana_corner
    use flow_data
    implicit none

    integer i, j, k, m, i1, ierr
    real(kind=OCFD_REAL_KIND)::   p00, tmp
    real(kind=OCFD_REAL_KIND), dimension(:), allocatable::  sx1, sy1, us1, Ts1, ps1, us0, Ts0, ps0

    p00 = 1.d0/(Para%gamma*Para%Ma*Para%Ma)
    allocate (sx1(nx), sy1(nx), us1(nx_global), Ts1(nx_global), ps1(nx_global), us0(nx_global), Ts0(nx_global), ps0(nx_global))

    us1 = 0.d0
    Ts1 = 0.d0
    ps1 = 0.d0
    us0 = 0.d0
    Ts0 = 0.d0
    ps0 = 0.d0

    if (npy .eq. 0) then

       if (npx .ne. npx0 - 1) then
          do i = 1, nx
             sx1(i) = Axx(i + 1, 1, 1) - Axx(i, 1, 1)
             sy1(i) = Ayy(i + 1, 1, 1) - Ayy(i, 1, 1)
          end do
       else
          do i = 1, nx - 1
             sx1(i) = Axx(i + 1, 1, 1) - Axx(i, 1, 1)
             sy1(i) = Ayy(i + 1, 1, 1) - Ayy(i, 1, 1)
          end do
          sx1(nx) = Axx(nx, 1, 1) - Axx(nx - 1, 1, 1)
          sy1(nx) = Ayy(nx, 1, 1) - Ayy(nx - 1, 1, 1)
       end if

       do i = 1, nx
          tmp = 1.d0/sqrt(sx1(i)**2 + sy1(i)**2)
          sx1(i) = sx1(i)*tmp
          sy1(i) = sy1(i)*tmp
       end do
!-----------------------
       do i = 1, nx
          i1 = i_offset(npx) + i - 1
          do k = 1, nz
             us1(i1) = us1(i1) + u(i, 2, k)*sx1(i) + v(i, 2, k)*sy1(i)
             Ts1(i1) = Ts1(i1) + T(i, 2, k)
             ps1(i1) = ps1(i1) + p00*d(i, 1, k)*T(i, 1, k)
          end do
       end do

    end if

    call MPI_REDUCE(us1, us0, nx_global, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(Ts1, Ts0, nx_global, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(ps1, ps0, nx_global, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    tmp = 1.d0/nz_global
    us0 = us0*tmp
    Ts0 = Ts0*tmp
    ps0 = ps0*tmp

    if (my_id .eq. 0) then
       open (101, file="WallUtp-log.dat", form="unformatted", position="append")
       write (101) Istep, tt
       write (101) us0, Ts0, ps0
       close (101)
    end if

    deallocate (sx1, sy1, us1, Ts1, ps1, us0, Ts0, ps0)

 end
!--------------------------------------------------------

