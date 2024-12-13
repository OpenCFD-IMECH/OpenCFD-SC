
!-----------------------------------------------------------------------------------
! 3 steps 3rd order TVD type Runge-Kutta method by Jiang & Shu
  subroutine OCFD_time_adv_RK3(KRK)
     use flow_data
     implicit none
     integer KRK, m, i, j, k
     real(kind=OCFD_REAL_KIND):: Ralfa(3), Rbeta(3)

     Ralfa(1) = 1.d0; Rbeta(1) = 1.d0
     Ralfa(2) = 3.d0/4.d0; Rbeta(2) = 1.d0/4.d0
     Ralfa(3) = 1.d0/3.d0; Rbeta(3) = 2.d0/3.d0

     if (KRK .eq. 1) then
        do m = 1, 5
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
!       f(i,j,m)=Ralfa(KRK)*fn(i,j,m) +dt*du(i,j,m)*Rbeta(KRK)
                 f(i, j, k, m) = fn(i, j, k, m) + Para%dt*du(i, j, k, m)
              end do
              end do
           end do
        end do
     else
        do m = 1, 5
           do k = 1, nz
              do j = 1, ny
              do i = 1, nx
                 f(i, j, k, m) = Ralfa(KRK)*fn(i, j, k, m) + Rbeta(KRK)*(f(i, j, k, m) + Para%dt*du(i, j, k, m))
              end do
              end do
           end do
        end do
     end if
  end

!----------------------------------------------

  subroutine comput_duvwT
     Use flow_data
     implicit none
     integer:: Num_NegT, i, j, k
     Num_NegT = 0
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        d(i, j, k) = f(i, j, k, 1)
        u(i, j, k) = f(i, j, k, 2)/f(i, j, k, 1)
        v(i, j, k) = f(i, j, k, 3)/f(i, j, k, 1)
        w(i, j, k) = f(i, j, k, 4)/f(i, j, k, 1)
        T(i,j,k)=(f(i,j,k,5) -(f(i,j,k,2)*u(i,j,k) +f(i,j,k,3)*v(i,j,k) +f(i,j,k,4)*w(i,j,k))*0.5d0 )/(f(i,j,k,1)*Cv)
!------------------------------------------------------------------
!  if T<=0, Error ! computation will stop (very useful message for debugging at overflow case)-------
        if (T(i, j, k) <= 0) then
           call handle_NegativeT(i, j, k, Num_NegT)
        end if
!----------------------------------
     end do
     end do
     end do

  end

  subroutine show_flow_msg(wall_time)
     use flow_data
     implicit none
     real*8:: wall_time, wtmp, E0(3), E1(3), p00, tmp0, tmp1
     integer:: i, j, k, ierr
     wtmp = wall_time
     wall_time = MPI_wtime()
     p00 = 1.d0/(Para%gamma*Para%Ma*Para%Ma)
     tmp1 = 1.d0 - Para%gamma

     E1(:) = 0.d0
     E0(:) = 0.d0
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        tmp0 = 1.d0/Ajac(i, j, k)
        E1(1) = E1(1) + f(i, j, k, 5)*tmp0                ! Total Energy
        E1(2) = E1(2) + d(i, j, k)*(u(i, j, k)*u(i, j, k) + v(i, j, k)*v(i, j, k) + w(i, j, k)*w(i, j, k))*0.5d0*tmp0  !Kinetic Energy
        E1(3) = E1(3) + p00*d(i, j, k)**tmp1*T(i, j, k)*tmp0                !Entropy  p/rho**gamma
     end do
     end do
     end do

     call MPI_REDUCE(E1, E0, 3, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     E0 = E0/(nx_global*ny_global*nz_global)

     if (my_id .eq. 0) then
        print *, 'Istep=', Istep, 'tt=', tt
        print *, 'CPU time per', Para%Istep_show, 'step is:    ', wall_time - wtmp
        write (*, "(3A25)") "Averaged Total-energy E","Kinetic-energy K","Entropy S"
        write (*, "(3E25.15)") E0(1), E0(2), E0(3)

        open (66, file='opencfd.log', position='append')
        write (66, *) 'Istep=', Istep, 'tt=', tt, 'CPU time is', wall_time - wtmp
        write (66, "(3A25)") "Averaged Total-energy E","Kinetic-energy K","Entropy S"
        write (66, "(3E25.15)") E0(1), E0(2), E0(3)
        close (66)
     end if

  end

