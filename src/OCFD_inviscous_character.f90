!c inviscous term for Jacobian- transformed N-S equation ---------------------
! Characteristic flux
  subroutine du_inviscous_Character
     use flow_data
     use inviscous_data
     implicit none

     integer:: i, j, k, m, Ka1, Kb1, Ka2, Kb2, set_hybrid_scheme
     real(kind=OCFD_REAL_KIND):: hx_1, hy_1, hz_1, Rhb0
!c-------------------------------------------------
     Ka1 = Scheme%Ka1
     Kb1 = Scheme%Kb1
     Ka2 = Scheme%Ka2
     Kb2 = Scheme%Kb2
     hx_1 = 1.d0/hx
     hy_1 = 1.d0/hy
     hz_1 = 1.d0/hz

     do k = 1 - LAP, nz + LAP
        do j = 1 - LAP, ny + LAP
        do i = 1 - LAP, nx + LAP
           cc(i, j, k) = sqrt(T(i, j, k))/Para%Ma      ! speed of sound
        end do
        end do
     end do

!c---------x direction-------------------------
     do k = 1, nz
     do j = 1, ny

        if (Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
           do i = 0, nx
              Rhb0 = max(Rhybrid(i, j, k), Rhybrid(i + 1, j, k))
              Scm_Hbx(i) = set_hybrid_scheme(Rhb0)
           end do
        end if

        do i = 1 - LAP, nx + LAP
           c1x(i) = cc(i, j, k)
           d1x(i) = d(i, j, k)
           u1x(i) = u(i, j, k)
           v1x(i) = v(i, j, k)
           w1x(i) = w(i, j, k)
           A1x(i) = Akx1(i, j, k)
           A2x(i) = Aky1(i, j, k)
           A3x(i) = Akz1(i, j, k)
        end do

        if (Para%Flux_Splitting .eq. OCFD_Split_SW) then
           call split_Stager_Warming(nx, LAP, c1x, d1x, u1x, v1x, w1x, A1x, A2x, A3x, fpx, fmx)
        else
           call split_Local_LaxFriedrichs(nx, LAP, c1x, d1x, u1x, v1x, w1x, A1x, A2x, A3x, fpx, fmx)
        end if

        call flux_charteric(nx, LAP, Ka1, Kb1, Ka2, Kb2, &
                            hhx, fpx, fmx, A1x, A2x, A3x, &
                            c1x, d1x, u1x, v1x, w1x, Scheme%Bound_index(1, 1), Scm_Hbx, Scheme%Scheme_boundary(1:2))
        do m = 1, 5
        do i = 1, nx
           du(i, j, k, m) = (hhx(i, m) - hhx(i - 1, m))*hx_1
        end do
        end do
     end do
     end do

!c-------y direction ---------------------------
     do k = 1, nz
     do i = 1, nx

        if (Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
           do j = 0, ny
              Rhb0 = max(Rhybrid(i, j, k), Rhybrid(i, j + 1, k))
              Scm_Hby(j) = set_hybrid_scheme(Rhb0)
           end do
        end if

        do j = 1 - LAP, ny + LAP
           c1y(j) = cc(i, j, k)
           d1y(j) = d(i, j, k)
           u1y(j) = u(i, j, k)
           v1y(j) = v(i, j, k)
           w1y(j) = w(i, j, k)
           A1y(j) = Aix1(i, j, k)
           A2y(j) = Aiy1(i, j, k)
           A3y(j) = Aiz1(i, j, k)
        end do

        if (Para%Flux_Splitting .eq. OCFD_Split_SW) then
           call split_Stager_Warming(ny, LAP, c1y, d1y, u1y, v1y, w1y, A1y, A2y, A3y, fpy, fmy)
        else
           call split_Local_LaxFriedrichs(ny, LAP, c1y, d1y, u1y, v1y, w1y, A1y, A2y, A3y, fpy, fmy)
        end if

        call flux_charteric(ny, LAP, Ka1, Kb1, Ka2, Kb2, &
                            hhy, fpy, fmy, A1y, A2y, A3y, &
                            c1y, d1y, u1y, v1y, w1y, Scheme%Bound_index(1, 2), Scm_Hby, Scheme%Scheme_boundary(3:4))

        do m = 1, 5
           do j = 1, ny
              du(i, j, k, m) = du(i, j, k, m) + (hhy(j, m) - hhy(j - 1, m))*hy_1
           end do
        end do
     end do
     end do
!c-------z direction ---------------------------
     do j = 1, ny
     do i = 1, nx

        if (Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
           do k = 0, nz
              Rhb0 = max(Rhybrid(i, j, k), Rhybrid(i, j, k + 1))
              Scm_Hbz(k) = set_hybrid_scheme(Rhb0)
           end do
        end if

        do k = 1 - LAP, nz + LAP
           c1z(k) = cc(i, j, k)
           d1z(k) = d(i, j, k)
           u1z(k) = u(i, j, k)
           v1z(k) = v(i, j, k)
           w1z(k) = w(i, j, k)
           A1z(k) = Asx1(i, j, k)
           A2z(k) = Asy1(i, j, k)
           A3z(k) = Asz1(i, j, k)
        end do

        if (Para%Flux_Splitting .eq. OCFD_Split_SW) then
           call split_Stager_Warming(nz, LAP, c1z, d1z, u1z, v1z, w1z, A1z, A2z, A3z, fpz, fmz)
        else
           call split_Local_LaxFriedrichs(nz, LAP, c1z, d1z, u1z, v1z, w1z, A1z, A2z, A3z, fpz, fmz)
        end if

        call flux_charteric(nz, LAP, Ka1, Kb1, Ka2, Kb2, &
                            hhz, fpz, fmz, A1z, A2z, A3z, &
                            c1z, d1z, u1z, v1z, w1z, Scheme%Bound_index(1, 3), Scm_Hbz, Scheme%Scheme_boundary(5:6))

        do m = 1, 5
           do k = 1, nz
              du(i, j, k, m) = -(du(i, j, k, m) + (hhz(k, m) - hhz(k - 1, m))*hz_1)*Ajac(i, j, k)                ! df/dt=-du
           end do
        end do
     end do
     end do

     if (Para%IF_Mass_Force .eq. 1) then
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
           du(i, j, k, 2) = du(i, j, k, 2) + d(i, j, k)*Para%Mass_Force(1)
           du(i, j, k, 3) = du(i, j, k, 3) + d(i, j, k)*Para%Mass_Force(2)
           du(i, j, k, 4) = du(i, j, k, 4) + d(i, j, k)*Para%Mass_Force(3)
           du(i, j, k, 5) = du(i, j, k, 5) + d(i, j, k)*(u(i, j, k)*Para%Mass_Force(1) + v(i, j, k)*Para%Mass_Force(2) &
                                                         + w(i, j, k)*Para%Mass_Force(3))
        end do
        end do
        end do
     end if

  end

!c-----------------------------------
  subroutine flux_charteric(nn, LAP1, Ka1, Kb1, Ka2, Kb2, hh, fp, fm, AA1, AA2, AA3, &
                            cc, dd, uu, vv, ww, Bound_index, Scm_Hb, Scheme_boundary)
     use flow_para
     implicit none

     integer nn, LAP1, Ka1, Kb1, Ka2, Kb2, i, m, mm, Bound_index(2), ia, ib, Scheme_boundary(2), Ka1h, Kb1h, Ka2h, Kb2h
     real(kind=OCFD_REAL_KIND), dimension(1 - LAP1:nn + LAP1):: cc, dd, uu, vv, ww, AA1, AA2, AA3
     real(kind=OCFD_REAL_KIND), dimension(1 - LAP1:nn + LAP1, 5):: fp, fm
     real(kind=OCFD_REAL_KIND):: hh(0:nn, 5)
     real(kind=OCFD_REAL_KIND):: fpc(Ka1:Kb1, 5), fmc(Ka2:Kb2, 5), h1(5), h2(5), h0(5)
     real(kind=OCFD_REAL_KIND):: S(5, 5), S1(5, 5)
     real(kind=OCFD_REAL_KIND):: d1, u1, v1, w1, c1, V2, un, ul, um
     real(kind=OCFD_REAL_KIND):: a1, a2, a3, ss, n1, n2, n3, l1, l2, l3, m1, m2, m3, KK1, KK, X1, H
     integer:: Scm_Hb(0:nn)        ! Hybrid scheme index (1 linear scheme ;  2 shock-capturing scheme)
     logical If_Char               ! Characteristic
!------------------------------------------------------
     if (Bound_index(1) .eq. 0 .or. Scheme_boundary(1) .ne. 0) then    ! 1 WENO5+ Ghost Cell ;  -1 Inner scheme + Ghost Cell
        ia = 0
     else
        ia = 1                  ! boundary scheme
     end if
     if (Bound_index(2) .eq. 0 .or. Scheme_boundary(2) .ne. 0) then
        ib = nn
     else
        ib = nn - 1
     end if
!-------------------------------------------------------

     do i = ia, ib

        if (Scm_Hb(i) == 1) then    ! linear scheme (DO NOT use character construction)
           fpc(Ka1:Kb1, :) = fp(i + Ka1:i + Kb1, :)
           fmc(Ka2:Kb2, :) = fm(i + Ka2:i + Kb2, :)
        else
           if (Scm_Hb(i) == 0) then    !  other nonlinear scheme (not hybrid scheme)
              Ka1h = Ka1; Kb1h = Kb1; Ka2h = Ka2; Kb2h = Kb2
           else if (Scm_Hb(i) == 2) then    ! WENO7 scheme
              Ka1h = Scheme%Ka1_H2; Kb1h = Scheme%Kb1_H2; Ka2h = Scheme%Ka2_H2; Kb2h = Scheme%Kb2_H2
           else                       ! WENO5 scheme
              Ka1h = Scheme%Ka1_H3; Kb1h = Scheme%Kb1_H3; Ka2h = Scheme%Ka2_H3; Kb2h = Scheme%Kb2_H3
           end if
           ! Shock capturing scheme (Using character construction)
           ! Transform variables to character space
           u1 = (uu(i) + uu(i + 1))*0.5d0
           v1 = (vv(i) + vv(i + 1))*0.5d0
           w1 = (ww(i) + ww(i + 1))*0.5d0
           c1 = (cc(i) + cc(i + 1))*0.5d0
           a1 = (AA1(i) + AA1(i + 1))*0.5d0
           a2 = (AA2(i) + AA2(i + 1))*0.5d0
           a3 = (AA3(i) + AA3(i + 1))*0.5d0

!  A=S(-1)*LAMDA*S   = R * LAMDA * L
!  See Katate Masatsuka's Book:  "I do like CFD, Vol. 1" page 77-78

           ss = sqrt(a1*a1 + a2*a2 + a3*a3)
           n1 = a1/ss; n2 = a2/ss; n3 = a3/ss

           if (abs(n3) <= abs(n2)) then
              ss = sqrt(n1*n1 + n2*n2)
              l1 = -n2/ss; l2 = n1/ss; l3 = 0.d0
           else
              ss = sqrt(n1*n1 + n3*n3)
              l1 = -n3/ss; l2 = 0.d0; l3 = n1/ss
           end if
           m1 = n2*l3 - n3*l2
           m2 = n3*l1 - n1*l3
           m3 = n1*l2 - n2*l1
           un = u1*n1 + v1*n2 + w1*n3
           ul = u1*l1 + v1*l2 + w1*l3
           um = u1*m1 + v1*m2 + w1*m3
           V2 = (u1*u1 + v1*v1 + w1*w1)*0.5d0
           KK = (para%gamma - 1.d0)/(c1*c1)
           KK1 = KK*0.5d0
           X1 = 1.d0/(2.d0*c1)
           H = V2 + 1.d0/KK
!==============S=L (Lift Characteristic Matrix)
           S(1, 1) = 1.d0 - KK*V2; S(1, 2) = KK*u1; S(1, 3) = KK*v1; S(1, 4) = KK*w1; S(1, 5) = -KK
           S(2, 1) = -ul; S(2, 2) = l1; S(2, 3) = l2; S(2, 4) = l3; S(2, 5) = 0.d0
           S(3, 1) = -um; S(3, 2) = m1; S(3, 3) = m2; S(3, 4) = m3; S(3, 5) = 0.d0
           S(4, 1) = KK1*V2 + X1*un; S(4, 2) = -X1*n1 - KK1*u1; S(4, 3) = -X1*n2 - KK1*v1; S(4, 4) = -X1*n3 - KK1*w1; S(4, 5) = KK1
           S(5, 1) = KK1*V2 - X1*un; S(5, 2) = X1*n1 - KK1*u1; S(5, 3) = X1*n2 - KK1*v1; S(5, 4) = X1*n3 - KK1*w1; S(5, 5) = KK1
!=======S1 = S^(-1)=R  (Right Characteristic Matrix)
           S1(1, 1) = 1.d0; S1(1, 2) = 0.d0; S1(1, 3) = 0.d0; S1(1, 4) = 1.d0; S1(1, 5) = 1.d0
           S1(2, 1) = u1; S1(2, 2) = l1; S1(2, 3) = m1; S1(2, 4) = u1 - c1*n1; S1(2, 5) = u1 + c1*n1
           S1(3, 1) = v1; S1(3, 2) = l2; S1(3, 3) = m2; S1(3, 4) = v1 - c1*n2; S1(3, 5) = v1 + c1*n2
           S1(4, 1) = w1; S1(4, 2) = l3; S1(4, 3) = m3; S1(4, 4) = w1 - c1*n3; S1(4, 5) = w1 + c1*n3
           S1(5, 1) = V2; S1(5, 2) = ul; S1(5, 3) = um; S1(5, 4) = H - c1*un; S1(5, 5) = H + c1*un

!         call test_S(S,S1)

! V=SU      V=S*F  ！ transform into character space
!       do mm=Ka1,Kb1
           do mm = Ka1h, Kb1h
              m = i + mm
              fpc(mm, 1) = S(1, 1)*fp(m, 1) + S(1, 2)*fp(m, 2) + S(1, 3)*fp(m, 3) + S(1, 4)*fp(m, 4) + S(1, 5)*fp(m, 5)
              fpc(mm, 2) = S(2, 1)*fp(m, 1) + S(2, 2)*fp(m, 2) + S(2, 3)*fp(m, 3) + S(2, 4)*fp(m, 4) + S(2, 5)*fp(m, 5)
              fpc(mm, 3) = S(3, 1)*fp(m, 1) + S(3, 2)*fp(m, 2) + S(3, 3)*fp(m, 3) + S(3, 4)*fp(m, 4) + S(3, 5)*fp(m, 5)
              fpc(mm, 4) = S(4, 1)*fp(m, 1) + S(4, 2)*fp(m, 2) + S(4, 3)*fp(m, 3) + S(4, 4)*fp(m, 4) + S(4, 5)*fp(m, 5)
              fpc(mm, 5) = S(5, 1)*fp(m, 1) + S(5, 2)*fp(m, 2) + S(5, 3)*fp(m, 3) + S(5, 4)*fp(m, 4) + S(5, 5)*fp(m, 5)
           end do

!       do mm=Ka2,Kb2
           do mm = Ka2h, Kb2h
              m = i + mm
              fmc(mm, 1) = S(1, 1)*fm(m, 1) + S(1, 2)*fm(m, 2) + S(1, 3)*fm(m, 3) + S(1, 4)*fm(m, 4) + S(1, 5)*fm(m, 5)
              fmc(mm, 2) = S(2, 1)*fm(m, 1) + S(2, 2)*fm(m, 2) + S(2, 3)*fm(m, 3) + S(2, 4)*fm(m, 4) + S(2, 5)*fm(m, 5)
              fmc(mm, 3) = S(3, 1)*fm(m, 1) + S(3, 2)*fm(m, 2) + S(3, 3)*fm(m, 3) + S(3, 4)*fm(m, 4) + S(3, 5)*fm(m, 5)
              fmc(mm, 4) = S(4, 1)*fm(m, 1) + S(4, 2)*fm(m, 2) + S(4, 3)*fm(m, 3) + S(4, 4)*fm(m, 4) + S(4, 5)*fm(m, 5)
              fmc(mm, 5) = S(5, 1)*fm(m, 1) + S(5, 2)*fm(m, 2) + S(5, 3)*fm(m, 3) + S(5, 4)*fm(m, 4) + S(5, 5)*fm(m, 5)
           end do
        end if

        do m = 1, 5
           ! scheme (including boundary scheme)
           call flux_x1(Ka1, Kb1, fpc(Ka1, m), h1(m), Scheme%Scheme_Invis, Bound_index, i, nn, Scm_Hb(i), Scheme_boundary)
           call flux_x2(Ka2, Kb2, fmc(Ka2, m), h2(m), Scheme%Scheme_Invis, Bound_index, i, nn, Scm_Hb(i), Scheme_boundary)
        end do

        h0 = h1 + h2
        if (Scm_Hb(i) == 1) then  ! linear scheme (do not use character flux)
           hh(i, :) = h0(:)
        else         ! Transform to phyics space
           hh(i, 1) = S1(1, 1)*h0(1) + S1(1, 2)*h0(2) + S1(1, 3)*h0(3) + S1(1, 4)*h0(4) + S1(1, 5)*h0(5)
           hh(i, 2) = S1(2, 1)*h0(1) + S1(2, 2)*h0(2) + S1(2, 3)*h0(3) + S1(2, 4)*h0(4) + S1(2, 5)*h0(5)
           hh(i, 3) = S1(3, 1)*h0(1) + S1(3, 2)*h0(2) + S1(3, 3)*h0(3) + S1(3, 4)*h0(4) + S1(3, 5)*h0(5)
           hh(i, 4) = S1(4, 1)*h0(1) + S1(4, 2)*h0(2) + S1(4, 3)*h0(3) + S1(4, 4)*h0(4) + S1(4, 5)*h0(5)
           hh(i, 5) = S1(5, 1)*h0(1) + S1(5, 2)*h0(2) + S1(5, 3)*h0(3) + S1(5, 4)*h0(4) + S1(5, 5)*h0(5)
        end if
     end do

!--------non Reflaction boundary (i==0,  i==nn) --------------
! 实现无反射边界条件， 在i=0处 hh+(0)=hh+ (1);  是的df+/dx =0  (i=0处）

     if (Bound_index(1) .eq. 1) then
        if (Scheme_boundary(1) .eq. 0) then  ! Default (WENO5-type)
           do m = 1, 5
              hh(0, m) = (2.d0*fp(1, m) + 5.d0*fp(2, m) - fp(3, m) + 2.d0*fm(3, m) - 7.d0*fm(2, m) + 11.d0*fm(1, m))/6.d0
           end do
        end if
     end if

     if (Bound_index(2) .eq. 1) then
        if (Scheme_boundary(2) .eq. 0) then  ! Default (WENO5-type)
        do m = 1, 5
  hh(nn, m) = (2.d0*fp(nn - 2, m) - 7.d0*fp(nn - 1, m) + 11.d0*fp(nn, m) + 2.d0*fm(nn, m) + 5.d0*fm(nn - 1, m) - fm(nn - 2, m))/6.d0
        end do
        end if
     end if
  end

!-----------------------------------------------------------------------------
  subroutine flux_x1(Ka, Kb, ff, hh, Num_Scheme, Bound_index, i, nn, Ihybrid, Scheme_boundary)
     use flow_para
     implicit none
     integer Num_Scheme, i, nn, Ka, Kb, Bound_index(2), Ihybrid, Scheme_boundary(2)
     real(kind=OCFD_REAL_KIND)::  ff(Ka:Kb), hh

     if ((Bound_index(1) .eq. 0 .or. i > -Ka .or. Scheme_boundary(1) .eq. -1) &
         .and. (Bound_index(2) .eq. 0 .or. i <= nn - Kb .or. Scheme_boundary(2) .eq. -1)) then
!     inner point
        if (Num_Scheme .eq. OCFD_Scheme_OMP6) then
           call hh_OMP6P(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_WENO7) then
           call hh_WENO7P(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_WENO5) then
           call hh_WENO5P(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_UD7L) then
           call hh_UD7L_P(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_Hybrid) then
           call hh_HybridP(Ka, Kb, ff, hh, Ihybrid)
        else if (Num_Scheme .eq. OCFD_Scheme_USER) then
           call hh_Scheme_USER_P(Ka, Kb, ff, hh)
        end if
        !  Boundary scheme
     else if (i <= -Ka .and. Bound_index(1) .eq. 1 .and. Scheme_boundary(1) .ne. -1) then
! Lift boundary point
        if (Scheme_boundary(1) .eq. 0) then  ! Default (WENO5-type)
           if (i == 1) then
              hh = (2.d0*ff(0) + 5.d0*ff(1) - ff(2))/6.d0               ! 3rd one-side
           else if (i == 2) then
              call hh_weno5P_boundary(Ka, Kb, ff, hh, DEL_LIFT)
           else
              call hh_weno5P(Ka, Kb, ff, hh)      ! full WENO5 scheme
           end if
        else if (Scheme_boundary(1) .eq. 1) then  ! WENO5 + Ghost Cell
           call hh_weno5P(Ka, Kb, ff, hh)      ! full WENO5 scheme
        end if

     else if (i > nn - Kb .and. Bound_index(2) .eq. 1 .and. Scheme_boundary(2) .ne. -1) then
! Right boundary point
        if (Scheme_boundary(2) .eq. 0) then  ! Default (WENO5-type)
           if (i == nn - 1) then
              call hh_weno5P_boundary(Ka, Kb, ff, hh, DEL_RIGHT)
           else
              call hh_WENO5P(Ka, Kb, ff, hh)
           end if
        else if (Scheme_boundary(2) .eq. 1) then  !  WENO5+ Ghost Cell
           call hh_WENO5P(Ka, Kb, ff, hh)
        end if

     end if

  end

!-----------------------------------------------------------------------------
  subroutine flux_x2(Ka, Kb, ff, hh, Num_Scheme, Bound_index, i, nn, Ihybrid, Scheme_boundary)
     use flow_para
     implicit none
     integer Num_Scheme, i, nn, Ka, Kb, Bound_index(2), Ihybrid, Scheme_boundary(2)
     real(kind=OCFD_REAL_KIND)::  ff(Ka:Kb), hh

     if ((Bound_index(1) .eq. 0 .or. i > -Ka .or. Scheme_boundary(1) .eq. -1) &
         .and. (Bound_index(2) .eq. 0 .or. i <= nn - Kb .or. Scheme_boundary(2) .eq. -1)) then
!     inner point
        if (Num_Scheme .eq. OCFD_Scheme_OMP6) then
           call hh_OMP6M(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_WENO7) then
           call hh_WENO7M(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_WENO5) then
           call hh_WENO5M(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_UD7L) then
           call hh_UD7L_M(Ka, Kb, ff, hh)
        else if (Num_Scheme .eq. OCFD_Scheme_Hybrid) then
           call hh_HybridM(Ka, Kb, ff, hh, Ihybrid)
        else if (Num_Scheme .eq. OCFD_Scheme_USER) then
           call hh_Scheme_USER_M(Ka, Kb, ff, hh)
        end if

     else if (i <= -Ka .and. Bound_index(1) .eq. 1 .and. Scheme_boundary(1) .ne. -1) then
        if (Scheme_boundary(1) .eq. 0) then  ! Default (WENO5-type)
           if (i == 1) then
              call hh_weno5M_boundary(Ka, Kb, ff, hh, DEL_LIFT)
           else
              call hh_WENO5M(Ka, Kb, ff, hh)
           end if
        else if (Scheme_boundary(1) .eq. 1) then  ! WENO5 + Ghost Cell
           call hh_WENO5M(Ka, Kb, ff, hh)
        end if

     else if (i > nn - Kb .and. Bound_index(2) .eq. 1 .and. Scheme_boundary(2) .ne. -1) then
        if (Scheme_boundary(2) .eq. 0) then  ! Default (WENO5-type)
           if (i == nn - 1) then
              hh = (2.d0*ff(1) + 5.d0*ff(0) - ff(-1))/6.d0
           else if (i == nn - 2) then
              call hh_weno5M_boundary(Ka, Kb, ff, hh, DEL_RIGHT)
           else
              call hh_WENO5M(Ka, Kb, ff, hh)
           end if
        else if (Scheme_boundary(2) .eq. 1) then  ! WENO5 + Ghost Cell
           call hh_WENO5M(Ka, Kb, ff, hh)
        end if

     end if

  end

!-------------------------------
! Test S and S1 (Characteristic Matrix)
  subroutine test_S(S, S1)
     implicit none
     integer:: i, j, m
     real*8:: S(5, 5), S1(5, 5), ST(5, 5)
     ST = 0.d0
     do j = 1, 5
     do i = 1, 5
        ST(i, j) = S(i, 1)*S1(1, j) + S(i, 2)*S1(2, j) + S(i, 3)*S1(3, j) + S(i, 4)*S1(4, j) + S(i, 5)*S1(5, j)
     end do
     end do
     print *, "---------test S ------------------"
     print *, " -----s=  "
     do i = 1, 5
        write (*, "(5F20.15)") S(i, :)
     end do
     print *, " -------S1=---"
     do i = 1, 5
        write (*, "(5F20.15)") S1(i, :)
     end do
     print *, "-------ST=----"
     do i = 1, 5
        write (*, "(5F20.15)") ST(i, :)
     end do
  end

