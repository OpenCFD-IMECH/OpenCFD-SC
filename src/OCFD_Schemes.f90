! Shock Capturing Schemes:  WENO5, WENO7,   OMP6
!   The same as that in OCFD2d_Scheme_1.f90,  but more efficiently;
!   与OCFD2d_Scheme_1.f90 中的代码算法与功能均相同， 且计算效率更高 （但不适用于特征分裂）
!   5th order WENO scheme   (for single scheme ib= 0, ie=nx )
! orient==1 for flux+ ;  -1 for flux-
     subroutine OCFD_weno5P(v, hh, nx, LAP, ib, ie)
        Use OCFD_constants
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
        real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
        real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

!      do i=0,nx
        do i = ib, ie
           S0 = 13.d0/12.d0*(v(i) - 2.d0*v(i + 1) + v(i + 2))**2 + 1.d0/4.d0*(3.d0*v(i) - 4.d0*v(i + 1) + v(i + 2))**2
           S1 = 13.d0/12.d0*(v(i - 1) - 2.d0*v(i) + v(i + 1))**2 + 1.d0/4.d0*(v(i - 1) - v(i + 1))**2
           S2 = 13.d0/12.d0*(v(i - 2) - 2.d0*v(i - 1) + v(i))**2 + 1.d0/4.d0*(v(i - 2) - 4.d0*v(i - 1) + 3.d0*v(i))**2
           a0 = C03/((ep + S0)**2)
           a1 = C13/((ep + S1)**2)
           a2 = C23/((ep + S2)**2)
           am = a0 + a1 + a2
           q03 = 1.d0/3.d0*v(i) + 5.d0/6.d0*v(i + 1) - 1.d0/6.d0*v(i + 2)
           q13 = -1.d0/6.d0*v(i - 1) + 5.d0/6.d0*v(i) + 1.d0/3.d0*v(i + 1)
           q23 = 1.d0/3.d0*v(i - 2) - 7.d0/6.d0*v(i - 1) + 11.d0/6.d0*v(i)
           hh(i) = (a0*q03 + a1*q13 + a2*q23)/am
        end do
     end

     subroutine OCFD_weno5M(v, hh, nx, LAP, ib, ie)
        Use OCFD_constants
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
        real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
        real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0
        ! for flux-
!      do i=1,nx+1
        do i = ib + 1, ie + 1
           S0 = 13.d0/12.d0*(v(i) - 2.d0*v(i - 1) + v(i - 2))**2 + 1.d0/4.d0*(3.d0*v(i) - 4.d0*v(i - 1) + v(i - 2))**2
           S1 = 13.d0/12.d0*(v(i + 1) - 2.d0*v(i) + v(i - 1))**2 + 1.d0/4.d0*(v(i + 1) - v(i - 1))**2
           S2 = 13.d0/12.d0*(v(i + 2) - 2.d0*v(i + 1) + v(i))**2 + 1.d0/4.d0*(v(i + 2) - 4.d0*v(i + 1) + 3.d0*v(i))**2
           a0 = C03/((ep + S0)**2)
           a1 = C13/((ep + S1)**2)
           a2 = C23/((ep + S2)**2)
           am = a0 + a1 + a2
           q03 = 1.d0/3.d0*v(i) + 5.d0/6.d0*v(i - 1) - 1.d0/6.d0*v(i - 2)
           q13 = -1.d0/6.d0*v(i + 1) + 5.d0/6.d0*v(i) + 1.d0/3.d0*v(i - 1)
           q23 = 1.d0/3.d0*v(i + 2) - 7.d0/6.d0*v(i + 1) + 11.d0/6.d0*v(i)
           hh(i - 1) = (a0*q03 + a1*q13 + a2*q23)/am
        end do

     end

!-----------------------------------------------------------------------------------------------------
!   7th order WENO-JS schemes

!----------------------------------------------------------------------------------
     subroutine OCFD_weno7P(v, hh, nx, LAP, ib, ie)
        Use OCFD_constants
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
        real(kind=OCFD_REAL_KIND)::  S0, S1, S2, S3, s10, s11, s12, s13, s20, s21, s22, s23, s30, s31, s32, s33, &
                                    a0, a1, a2, a3, am, q0, q1, q2, q3
        real(kind=OCFD_REAL_KIND), parameter:: &
           C0 = 1.d0/35.d0, C1 = 12.d0/35.d0, C2 = 18.d0/35.d0, C3 = 4.d0/35.d0, &
           a11 = -2.d0/6.d0, a12 = 9.d0/6.d0, a13 = -18.d0/6.d0, a14 = 11.d0/6.d0, &
           a21 = 1.d0/6.d0, a23 = 3.d0/6.d0, a24 = 2.d0/6.d0, &
           a31 = -2.d0/6.d0, a32 = -3.d0/6.d0, a34 = -1.d0/6.d0, &
           a41 = -11.d0/6.d0, a42 = 18.d0/6.d0, a43 = -9.d0/6.d0, a44 = 2.d0/6.d0, &
           b12 = 4.d0, b13 = -5.d0, b14 = 2.d0, b22 = -2.d0, &
           b41 = 2.d0, b42 = -5.d0, b43 = 4.d0, c12 = 3.d0, &
           d12 = 13.d0/12.d0, d13 = 1043.d0/960.d0, d14 = 1.d0/12.d0
        real(kind=OCFD_REAL_KIND), parameter:: &
           e11 = -3.d0/12.d0, e12 = 13.d0/12.d0, e13 = -23.d0/12.d0, e14 = 25.d0/12.d0, &
           e21 = 1.d0/12.d0, e22 = -5.d0/12.d0, e23 = 13.d0/12.d0, e24 = 3.d0/12.d0, &
           e31 = -1.d0/12.d0, e32 = 7.d0/12.d0, e33 = 7.d0/12.d0, e34 = -1.d0/12.d0, &
           e41 = 3.d0/12.d0, e42 = 13.d0/12.d0, e43 = -5.d0/12.d0, e44 = 1.d0/12.d0, &
           ep = 1.d-8    !! WENO-JS

!  do i=0,nx
        do i = ib, ie
! 7th order WENO scheme
! 1  阶导数
           S10 = a11*v(i - 3) + a12*v(i - 2) + a13*v(i - 1) + a14*v(i)
           S11 = a21*v(i - 2) - v(i - 1) + a23*v(i) + a24*v(i + 1)
           S12 = a31*v(i - 1) + a32*v(i) + v(i + 1) + a34*v(i + 2)
           S13 = a41*v(i) + a42*v(i + 1) + a43*v(i + 2) + a44*v(i + 3)
           ! 2 阶导数
           S20 = -v(i - 3) + b12*v(i - 2) + b13*v(i - 1) + b14*v(i)
           S21 = v(i - 1) + b22*v(i) + v(i + 1)
           S22 = v(i) + b22*v(i + 1) + v(i + 2)
           S23 = b41*v(i) + b42*v(i + 1) + b43*v(i + 2) - v(i + 3)
! 3 阶导数
           S30 = -v(i - 3) + c12*(v(i - 2) - v(i - 1)) + v(i)
           S31 = -v(i - 2) + c12*(v(i - 1) - v(i)) + v(i + 1)
           S32 = -v(i - 1) + c12*(v(i) - v(i + 1)) + v(i + 2)
           S33 = -v(i) + c12*(v(i + 1) - v(i + 2)) + v(i + 3)

           S0 = S10*S10 + d12*S20*S20 + d13*S30*S30 + d14*S10*S30
           S1 = S11*S11 + d12*S21*S21 + d13*S31*S31 + d14*S11*S31
           S2 = S12*S12 + d12*S22*S22 + d13*S32*S32 + d14*S12*S32
           S3 = S13*S13 + d12*S23*S23 + d13*S33*S33 + d14*S13*S33

!-------WENO J-S----------------------
           a0 = C0/((ep + S0)**2)
           a1 = C1/((ep + S1)**2)
           a2 = C2/((ep + S2)**2)
           a3 = C3/((ep + S3)**2)
!-----------------------------------------------
           am = a0 + a1 + a2 + a3

!  4阶差分格式的通量
           q0 = e11*v(i - 3) + e12*v(i - 2) + e13*v(i - 1) + e14*v(i)
           q1 = e21*v(i - 2) + e22*v(i - 1) + e23*v(i) + e24*v(i + 1)
           q2 = e31*v(i - 1) + e32*v(i) + e33*v(i + 1) + e34*v(i + 2)
           q3 = e41*v(i) + e42*v(i + 1) + e43*v(i + 2) + e44*v(i + 3)

!  由4个4阶差分格式组合成1个7阶差分格式
!     hj(i)=W0*q0+W1*q1+W2*q2+W3*q3
           hh(i) = (a0*q0 + a1*q1 + a2*q2 + a3*q3)/am
        end do
     end

!----------WENO7 for flux-
     subroutine OCFD_weno7M(v, hh, nx, LAP, ib, ie)
        Use OCFD_constants
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
        real(kind=OCFD_REAL_KIND)::  S0, S1, S2, S3, s10, s11, s12, s13, s20, s21, s22, s23, s30, s31, s32, s33, &
                                    a0, a1, a2, a3, am, q0, q1, q2, q3
        real(kind=OCFD_REAL_KIND), parameter:: &
           C0 = 1.d0/35.d0, C1 = 12.d0/35.d0, C2 = 18.d0/35.d0, C3 = 4.d0/35.d0, &
           a11 = -2.d0/6.d0, a12 = 9.d0/6.d0, a13 = -18.d0/6.d0, a14 = 11.d0/6.d0, &
           a21 = 1.d0/6.d0, a23 = 3.d0/6.d0, a24 = 2.d0/6.d0, &
           a31 = -2.d0/6.d0, a32 = -3.d0/6.d0, a34 = -1.d0/6.d0, &
           a41 = -11.d0/6.d0, a42 = 18.d0/6.d0, a43 = -9.d0/6.d0, a44 = 2.d0/6.d0, &
           b12 = 4.d0, b13 = -5.d0, b14 = 2.d0, b22 = -2.d0, &
           b41 = 2.d0, b42 = -5.d0, b43 = 4.d0, c12 = 3.d0, &
           d12 = 13.d0/12.d0, d13 = 1043.d0/960.d0, d14 = 1.d0/12.d0
        real(kind=OCFD_REAL_KIND), parameter:: &
           e11 = -3.d0/12.d0, e12 = 13.d0/12.d0, e13 = -23.d0/12.d0, e14 = 25.d0/12.d0, &
           e21 = 1.d0/12.d0, e22 = -5.d0/12.d0, e23 = 13.d0/12.d0, e24 = 3.d0/12.d0, &
           e31 = -1.d0/12.d0, e32 = 7.d0/12.d0, e33 = 7.d0/12.d0, e34 = -1.d0/12.d0, &
           e41 = 3.d0/12.d0, e42 = 13.d0/12.d0, e43 = -5.d0/12.d0, e44 = 1.d0/12.d0, &
           ep = 1.d-8    !! WENO-JS

!    do i= 1,nx+1
        do i = ib + 1, ie + 1
!      7th order WENO scheme
! 1  阶导数
           S10 = a11*v(i + 3) + a12*v(i + 2) + a13*v(i + 1) + a14*v(i)
           S11 = a21*v(i + 2) - v(i + 1) + a23*v(i) + a24*v(i - 1)
           S12 = a31*v(i + 1) + a32*v(i) + v(i - 1) + a34*v(i - 2)
           S13 = a41*v(i) + a42*v(i - 1) + a43*v(i - 2) + a44*v(i - 3)
! 2 阶导数
           S20 = -v(i + 3) + b12*v(i + 2) + b13*v(i + 1) + b14*v(i)
           S21 = v(i + 1) + b22*v(i) + v(i - 1)
           S22 = v(i) + b22*v(i - 1) + v(i - 2)
           S23 = b41*v(i) + b42*v(i - 1) + b43*v(i - 2) - v(i - 3)
! 3 阶导数
           S30 = -v(i + 3) + c12*(v(i + 2) - v(i + 1)) + v(i)
           S31 = -v(i + 2) + c12*(v(i + 1) - v(i)) + v(i - 1)
           S32 = -v(i + 1) + c12*(v(i) - v(i - 1)) + v(i - 2)
           S33 = -v(i) + c12*(v(i - 1) - v(i - 2)) + v(i - 3)

           S0 = S10*S10 + d12*S20*S20 + d13*S30*S30 + d14*S10*S30
           S1 = S11*S11 + d12*S21*S21 + d13*S31*S31 + d14*S11*S31
           S2 = S12*S12 + d12*S22*S22 + d13*S32*S32 + d14*S12*S32
           S3 = S13*S13 + d12*S23*S23 + d13*S33*S33 + d14*S13*S33

           a0 = C0/((ep + S0)**2)
           a1 = C1/((ep + S1)**2)
           a2 = C2/((ep + S2)**2)
           a3 = C3/((ep + S3)**2)

!-----------------------------------------------

           am = a0 + a1 + a2 + a3

!  4阶差分格式的通量
           q0 = e11*v(i + 3) + e12*v(i + 2) + e13*v(i + 1) + e14*v(i)
           q1 = e21*v(i + 2) + e22*v(i + 1) + e23*v(i) + e24*v(i - 1)
           q2 = e31*v(i + 1) + e32*v(i) + e33*v(i - 1) + e34*v(i - 2)
           q3 = e41*v(i) + e42*v(i - 1) + e43*v(i - 2) + e44*v(i - 3)

!  由4个4阶差分格式组合成1个7阶差分格式
           hh(i - 1) = (a0*q0 + a1*q1 + a2*q2 + a3*q3)/am
        end do

     end

!---------------------------------------------

!---------------------------------------------
! Optimized 6th order Monotonicity-Preserving Schemes
! Scheme by Leng Yan & Li Xinliang  see: Int. J. Numer. Meth. Fluids 2013; 73:560–577
!================================================================================
     subroutine OCFD_OMP6P(v, hh, nx, LAP, ib, ie)
        Use OCFD_constants
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
        real(kind=OCFD_REAL_KIND)::mid_nf
        real(kind=OCFD_REAL_KIND)::minmod2, minmod4
        real(kind=OCFD_REAL_KIND)::d1, d2, d3, ul, md, lc, mp, fmax, fmin
        real(kind=OCFD_REAL_KIND), parameter ::kappa = 4.d0, ep = 1.e-10
        real(kind=OCFD_REAL_KIND), parameter::  OMP_m = 0.015d0, OMP_n = 0.d0   !defaut (robust) ; !!! low viscous:  OMP_m=0.001d0 ; OMP_n=0.d0
        real(kind=OCFD_REAL_KIND), parameter:: OMP_a = 0.5d0*(OMP_m + OMP_n), OMP_b = 0.5d0*(OMP_m - OMP_n)
        real(kind=OCFD_REAL_KIND), parameter:: OMP_a1 = OMP_a, OMP_a2 = 1.d0/60.d0 - OMP_b - 6.d0*OMP_a, &
                                               OMP_a3 = -2.d0/15.d0 + 6.d0*OMP_b + 15.d0*OMP_a, &
                               OMP_a4 = 37.d0/60.d0 - 15.d0*OMP_b - 20.d0*OMP_a, OMP_a5 = 37.d0/60.d0 + 20.d0*OMP_b + 15.d0*OMP_a, &
                          OMP_a6 = -2.d0/15.d0 - 6.d0*OMP_a - 15.d0*OMP_b, OMP_a7 = 1.d0/60.d0 + OMP_a + 6.d0*OMP_b, OMP_a8 = -OMP_b

!        do i=0,nx
        do i = ib, ie
           mid_nf=OMP_a1*v(i+4)+OMP_a2*v(i+3)+OMP_a3*v(i+2)+OMP_a4*v(i+1)+OMP_a5*v(i)+OMP_a6*v(i-1)+OMP_a7*v(i-2)+OMP_a8*v(i-3)
           mp = v(i) + minmod2((v(i + 1) - v(i)), kappa*(v(i) - v(i - 1)))
           if ((mid_nf - v(i))*(mid_nf - mp) .ge. ep) then
              d1 = v(i - 2) + v(i) - 2.d0*v(i - 1)
              d2 = v(i - 1) + v(i + 1) - 2.d0*v(i)
              d3 = v(i) + v(i + 2) - 2.d0*v(i + 1)
              ul = v(i) + kappa*(v(i) - v(i - 1))
              md = 0.5d0*(v(i) + v(i + 1)) - 0.5d0*minmod4(4.d0*d2 - d3, 4.d0*d3 - d2, d2, d3)
              lc = v(i) + 0.5d0*(v(i) - v(i - 1)) + kappa*minmod4(4.d0*d1 - d2, 4.d0*d2 - d1, d2, d1)/3.d0
              fmin = max(min(v(i), v(i + 1), md), min(v(i), ul, lc))
              fmax = min(max(v(i), v(i + 1), md), max(v(i), ul, lc))
              mid_nf = mid_nf + minmod2(fmax - mid_nf, fmin - mid_nf)
           end if
           hh(i) = mid_nf
        end do
     end

     subroutine OCFD_OMP6M(v, hh, nx, LAP, ib, ie)
        Use OCFD_constants
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
        real(kind=OCFD_REAL_KIND)::mid_nf
        real(kind=OCFD_REAL_KIND)::minmod2, minmod4
        real(kind=OCFD_REAL_KIND)::d1, d2, d3, ul, md, lc, mp, fmax, fmin
        real(kind=OCFD_REAL_KIND), parameter ::kappa = 4.d0, ep = 1.e-10
        real(kind=OCFD_REAL_KIND), parameter::  OMP_m = 0.015d0, OMP_n = 0.d0   !defaut (robust) ; !!! low viscous:  OMP_m=0.001d0 ; OMP_n=0.d0
        real(kind=OCFD_REAL_KIND), parameter:: OMP_a = 0.5d0*(OMP_m + OMP_n), OMP_b = 0.5d0*(OMP_m - OMP_n)
        real(kind=OCFD_REAL_KIND), parameter:: OMP_a1 = OMP_a, OMP_a2 = 1.d0/60.d0 - OMP_b - 6.d0*OMP_a, &
                                               OMP_a3 = -2.d0/15.d0 + 6.d0*OMP_b + 15.d0*OMP_a, &
                               OMP_a4 = 37.d0/60.d0 - 15.d0*OMP_b - 20.d0*OMP_a, OMP_a5 = 37.d0/60.d0 + 20.d0*OMP_b + 15.d0*OMP_a, &
                          OMP_a6 = -2.d0/15.d0 - 6.d0*OMP_a - 15.d0*OMP_b, OMP_a7 = 1.d0/60.d0 + OMP_a + 6.d0*OMP_b, OMP_a8 = -OMP_b

!          do i=0,nx
        do i = ib, ie
           mid_nf=OMP_a1*v(i-3)+OMP_a2*v(i-2)+OMP_a3*v(i-1)+OMP_a4*v(i) +OMP_a5*v(i+1)+OMP_a6*v(i+2)+OMP_a7*v(i+3)+OMP_a8*v(i+4)
           mp = v(i + 1) + minmod2((v(i) - v(i + 1)), kappa*(v(i + 1) - v(i + 2)))
           if ((mid_nf - v(i + 1))*(mid_nf - mp) .ge. ep) then
              d1 = v(i + 3) + v(i + 1) - 2.d0*v(i + 2)
              d2 = v(i + 2) + v(i) - 2.d0*v(i + 1)
              d3 = v(i + 1) + v(i - 1) - 2.d0*v(i)
              !----
              ul = v(i + 1) + kappa*(v(i + 1) - v(i + 2))
              md = 0.5d0*(v(i + 1) + v(i)) - 0.5d0*minmod4(4.d0*d2 - d3, 4.d0*d3 - d2, d2, d3)
              lc = v(i + 1) + 0.5d0*(v(i + 1) - v(i)) + kappa*minmod4(4.d0*d1 - d2, 4.d0*d2 - d1, d2, d1)/3.d0
              fmin = max(min(v(i + 1), v(i), md), min(v(i + 1), ul, lc))
              fmax = min(max(v(i + 1), v(i), md), max(v(i + 1), ul, lc))
              mid_nf = mid_nf + minmod2(fmax - mid_nf, fmin - mid_nf)
           end if
           hh(i) = mid_nf
        end do

     end

! ==========================================
     function minmod2(x1, x2)
        Use OCFD_precision
        implicit none
        real(kind=OCFD_REAL_KIND):: x1, x2, minmod2
        minmod2 = 0.5d0*(sign(1.d0, x1) + sign(1.d0, x2))*min(abs(x1), abs(x2))
     end
!=========================================================
     function minmod4(x1, x2, x3, x4)
        Use OCFD_precision
        implicit none
        real(kind=OCFD_REAL_KIND):: x1, x2, x3, x4, minmod4
        minmod4 = 0.5d0*(sign(1.d0, x1) + sign(1.d0, x2))
        minmod4 = minmod4*abs(0.5d0*(sign(1.d0, x1) + sign(1.d0, x3)))
        minmod4 = minmod4*abs(0.5d0*(sign(1.d0, x1) + sign(1.d0, x4)))
        minmod4 = minmod4*min(abs(x1), abs(x2), abs(x3), abs(x4))
     end

!-------7th order low-dissipation scheme (mixing 7th upwind and 8th center)
     subroutine OCFD_UD7L_P(v, hh, nx, LAP, ib, ie)
        use OCFD_constants
        use Scheme_para
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
!   do i=0,nx
        do i = ib, ie
           hh(i) = Scheme%UD7L(1)*v(i - 3) + Scheme%UD7L(2)*v(i - 2) + Scheme%UD7L(3)*v(i - 1) + Scheme%UD7L(4)*v(i) &
                   + Scheme%UD7L(5)*v(i + 1) + Scheme%UD7L(6)*v(i + 2) + Scheme%UD7L(7)*v(i + 3) + Scheme%UD7L(8)*v(i + 4)
        end do
     end

     subroutine OCFD_UD7L_M(v, hh, nx, LAP, ib, ie)
        use OCFD_constants
        use Scheme_para
        implicit none
        integer nx, LAP, i, ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)
!   do i=0,nx
        do i = ib, ie
           hh(i) = Scheme%UD7L(1)*v(i + 4) + Scheme%UD7L(2)*v(i + 3) + Scheme%UD7L(3)*v(i + 2) + Scheme%UD7L(4)*v(i + 1) &
                   + Scheme%UD7L(5)*v(i) + Scheme%UD7L(6)*v(i - 1) + Scheme%UD7L(7)*v(i - 2) + Scheme%UD7L(8)*v(i - 3)
        end do
     end

!-------Hybrid Scheme ------------( HYBRID with UD7L and WENO5)
     subroutine OCFD_HybridP(v, hh, nx, LAP, Scm_Hy, ib, ie)
        use OCFD_constants
        use Scheme_para
        implicit none
        integer nx, LAP, i, Scm_Hy(0:nx), ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)

!    do i=0,nx
        do i = ib, ie
        if (Scm_Hy(i) == 1) then
           hh(i) = Scheme%UD7L(1)*v(i - 3) + Scheme%UD7L(2)*v(i - 2) + Scheme%UD7L(3)*v(i - 1) + Scheme%UD7L(4)*v(i) &
                   + Scheme%UD7L(5)*v(i + 1) + Scheme%UD7L(6)*v(i + 2) + Scheme%UD7L(7)*v(i + 3) + Scheme%UD7L(8)*v(i + 4)
        else if (Scm_Hy(i) == 2) then   ! WENO 7P
           call hh_weno7P(-3, 3, v(i - 3), hh(i))     ! v(i-3) === v(i-3:i+3)
        else
           call hh_weno5P(-2, 2, v(i - 2), hh(i))
        end if
        end do
     end

     subroutine OCFD_HybridM(v, hh, nx, LAP, Scm_Hy, ib, ie)
        Use OCFD_constants
        use Scheme_para
        implicit none
        integer nx, LAP, i, Scm_Hy(0:nx), ib, ie
        real(kind=OCFD_REAL_KIND):: v(1 - LAP:nx + LAP), hh(0:nx)

        ! for flux-
!    do i=0,nx
        do i = ib, ie
        if (Scm_Hy(i) == 1) then
           hh(i) = Scheme%UD7L(1)*v(i + 4) + Scheme%UD7L(2)*v(i + 3) + Scheme%UD7L(3)*v(i + 2) + Scheme%UD7L(4)*v(i + 1) &
                   + Scheme%UD7L(5)*v(i) + Scheme%UD7L(6)*v(i - 1) + Scheme%UD7L(7)*v(i - 2) + Scheme%UD7L(8)*v(i - 3)
        else if (Scm_Hy(i) == 2) then
           call hh_weno7M(-2, 4, v(i - 2), hh(i))
        else
           call hh_weno5M(-1, 3, v(i - 1), hh(i))
        end if
        end do
     end
