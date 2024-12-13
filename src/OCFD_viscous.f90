!========================================================
! Viscous term   ! 3D Jacobian transformation
! Copyright by Li Xinliang
module viscous_data                 ! data used by inviscous term
   use OCFD_precision
   implicit none
   real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :, :, :)::  Ev1, Ev2, Ev3
   real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :, :)::  uk, ui, us, vk, vi, vs, wk, wi, ws, Tk, Ti, Ts
end

subroutine allocate_vicous_data
   Use flow_para
   Use viscous_data
   implicit none
   allocate (Ev1(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP, 4), &
             Ev2(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP, 4), &
             Ev3(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP, 4))
   allocate (uk(nx, ny, nz), ui(nx, ny, nz), us(nx, ny, nz), vk(nx, ny, nz), vi(nx, ny, nz), vs(nx, ny, nz), &
             wk(nx, ny, nz), wi(nx, ny, nz), ws(nx, ny, nz), Tk(nx, ny, nz), Ti(nx, ny, nz), Ts(nx, ny, nz))
end

subroutine deallocate_vicous_data
   Use viscous_data
   implicit none
   deallocate (Ev1, Ev2, Ev3, uk, ui, us, vk, vi, vs, wk, wi, ws, Tk, Ti, Ts)
end

!-------------Viscous term -----------------------------------------------------------------
subroutine du_viscous
   Use flow_data
   Use viscous_data
   implicit none

   real(kind=OCFD_REAL_KIND)::  div, ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz, Amu1, Amuk, &
                               s11, s12, s13, s22, s23, s33, E1, E2, E3
   real(kind=OCFD_REAL_KIND), parameter:: Prt = 0.9d0, tmp2_3 = 2.d0/3.d0
   integer:: i, j, k, m
!------------------------------------------------------------------------

   call comput_Amu         ! comput viscous coefficient (by using Sutherland eq.)

   call OCFD_dx0(u, uk, Scheme%Scheme_Vis)
   call OCFD_dx0(v, vk, Scheme%Scheme_Vis)
   call OCFD_dx0(w, wk, Scheme%Scheme_Vis)
   call OCFD_dx0(T, Tk, Scheme%Scheme_Vis)
   call OCFD_dy0(u, ui, Scheme%Scheme_Vis)
   call OCFD_dy0(v, vi, Scheme%Scheme_Vis)
   call OCFD_dy0(w, wi, Scheme%Scheme_Vis)
   call OCFD_dy0(T, Ti, Scheme%Scheme_Vis)
   call OCFD_dz0(u, us, Scheme%Scheme_Vis)
   call OCFD_dz0(v, vs, Scheme%Scheme_Vis)
   call OCFD_dz0(w, ws, Scheme%Scheme_Vis)
   call OCFD_dz0(T, Ts, Scheme%Scheme_Vis)

!-------------------------------------------------------------

   do k = 1, nz
   do j = 1, ny
   do i = 1, nx
      ux = uk(i, j, k)*Akx(i, j, k) + ui(i, j, k)*Aix(i, j, k) + us(i, j, k)*Asx(i, j, k)
      vx = vk(i, j, k)*Akx(i, j, k) + vi(i, j, k)*Aix(i, j, k) + vs(i, j, k)*Asx(i, j, k)
      wx = wk(i, j, k)*Akx(i, j, k) + wi(i, j, k)*Aix(i, j, k) + ws(i, j, k)*Asx(i, j, k)
      Tx = Tk(i, j, k)*Akx(i, j, k) + Ti(i, j, k)*Aix(i, j, k) + Ts(i, j, k)*Asx(i, j, k)

      uy = uk(i, j, k)*Aky(i, j, k) + ui(i, j, k)*Aiy(i, j, k) + us(i, j, k)*Asy(i, j, k)
      vy = vk(i, j, k)*Aky(i, j, k) + vi(i, j, k)*Aiy(i, j, k) + vs(i, j, k)*Asy(i, j, k)
      wy = wk(i, j, k)*Aky(i, j, k) + wi(i, j, k)*Aiy(i, j, k) + ws(i, j, k)*Asy(i, j, k)
      Ty = Tk(i, j, k)*Aky(i, j, k) + Ti(i, j, k)*Aiy(i, j, k) + Ts(i, j, k)*Asy(i, j, k)

      uz = uk(i, j, k)*Akz(i, j, k) + ui(i, j, k)*Aiz(i, j, k) + us(i, j, k)*Asz(i, j, k)
      vz = vk(i, j, k)*Akz(i, j, k) + vi(i, j, k)*Aiz(i, j, k) + vs(i, j, k)*Asz(i, j, k)
      wz = wk(i, j, k)*Akz(i, j, k) + wi(i, j, k)*Aiz(i, j, k) + ws(i, j, k)*Asz(i, j, k)
      Tz = Tk(i, j, k)*Akz(i, j, k) + Ti(i, j, k)*Aiz(i, j, k) + Ts(i, j, k)*Asz(i, j, k)
      div = ux + vy + wz

!        Amu1=Amu(i,j,k)+Amu_t(i,j,k)   ! Turbulent or LES model !
!        Amuk=Cp*(Amu(i,j,k)/Pr+Amu_t(i,j,k)/Prt)

      Amu1 = Amu(i, j, k)            !  In this version, turbulence model is not supported !
      Amuk = Cp*Amu(i, j, k)/Para%Pr

      s11 = (2.d0*ux - tmp2_3*div)*Amu1          ! tmp2_3=2.d0/3.d0
      s22 = (2.d0*vy - tmp2_3*div)*Amu1
      s33 = (2.d0*wz - tmp2_3*div)*Amu1

      s12 = (uy + vx)*Amu1
      s13 = (uz + wx)*Amu1
      s23 = (vz + wy)*Amu1

      E1 = u(i, j, k)*s11 + v(i, j, k)*s12 + w(i, j, k)*s13 + Amuk*Tx
      E2 = u(i, j, k)*s12 + v(i, j, k)*s22 + w(i, j, k)*s23 + Amuk*Ty
      E3 = u(i, j, k)*s13 + v(i, j, k)*s23 + w(i, j, k)*s33 + Amuk*Tz

      Ev1(i, j, k, 1) = Akx1(i, j, k)*s11 + Aky1(i, j, k)*s12 + Akz1(i, j, k)*s13
      Ev1(i, j, k, 2) = Akx1(i, j, k)*s12 + Aky1(i, j, k)*s22 + Akz1(i, j, k)*S23
      Ev1(i, j, k, 3) = Akx1(i, j, k)*s13 + Aky1(i, j, k)*s23 + Akz1(i, j, k)*s33
      Ev1(i, j, k, 4) = Akx1(i, j, k)*E1 + Aky1(i, j, k)*E2 + Akz1(i, j, k)*E3

      Ev2(i, j, k, 1) = Aix1(i, j, k)*s11 + Aiy1(i, j, k)*s12 + Aiz1(i, j, k)*s13
      Ev2(i, j, k, 2) = Aix1(i, j, k)*s12 + Aiy1(i, j, k)*s22 + Aiz1(i, j, k)*S23
      Ev2(i, j, k, 3) = Aix1(i, j, k)*s13 + Aiy1(i, j, k)*s23 + Aiz1(i, j, k)*s33
      Ev2(i, j, k, 4) = Aix1(i, j, k)*E1 + Aiy1(i, j, k)*E2 + Aiz1(i, j, k)*E3

      Ev3(i, j, k, 1) = Asx1(i, j, k)*s11 + Asy1(i, j, k)*s12 + Asz1(i, j, k)*s13
      Ev3(i, j, k, 2) = Asx1(i, j, k)*s12 + Asy1(i, j, k)*s22 + Asz1(i, j, k)*S23
      Ev3(i, j, k, 3) = Asx1(i, j, k)*s13 + Asy1(i, j, k)*s23 + Asz1(i, j, k)*s33
      Ev3(i, j, k, 4) = Asx1(i, j, k)*E1 + Asy1(i, j, k)*E2 + Asz1(i, j, k)*E3

   end do
   end do
   end do

!------  x,y,z direction
   do m = 1, 4
      call exchange_boundary_x(Ev1(1 - LAP, 1 - LAP, 1 - LAP, m))
      call exchange_boundary_y(Ev2(1 - LAP, 1 - LAP, 1 - LAP, m))
      call exchange_boundary_z(Ev3(1 - LAP, 1 - LAP, 1 - LAP, m))

      call OCFD_dx0(Ev1(1 - LAP, 1 - LAP, 1 - LAP, m), uk, Scheme%Scheme_Vis)
      call OCFD_dy0(Ev2(1 - LAP, 1 - LAP, 1 - LAP, m), ui, Scheme%Scheme_Vis)
      call OCFD_dz0(Ev3(1 - LAP, 1 - LAP, 1 - LAP, m), us, Scheme%Scheme_Vis)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         du(i, j, k, m + 1) = du(i, j, k, m + 1) + (uk(i, j, k) + ui(i, j, k) + us(i, j, k))*Ajac(i, j, k)
      end do
      end do
      end do
   end do

end

!----------------------------------------------------------
! viscous coefficient: Sutherland Eq.
subroutine comput_Amu
   use flow_data
   implicit none
   real(kind=OCFD_REAL_KIND):: Tsb, tmpR
   integer:: i, j, k
   Tsb = 110.4d0/Para%Ref_Amu_T0
   TmpR = 1.d0/Para%Re
   do k = 1, nz
   do j = 1, ny
   do i = 1, nx
      Amu(i, j, k) = TmpR*(1.d0 + Tsb)*sqrt(T(i, j, k)**3)/(Tsb + T(i, j, k))
   end do
   end do
   end do
end

!c========================================================
!c      subroutine add_symmetry_Ev(kv,Ev,nx,ny,nz,LAP)
!c        include 'OpenCFD.h'
!c      real Ev
!c      integer i,j,k,LAP,kv
!c        dimension Ev(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
!c        end
!c--------------------------------------------------------
