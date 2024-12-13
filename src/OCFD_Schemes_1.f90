! Shock Capturing Schemes:  WENO5, WENO7,   OMP6
! 实现的功能于OCFD_Scheme.f90 中的程序相同， 算法也相同， 只是形式有些差异，更易于用于特征通量分裂  （但效率有所降低）

!--------hh :  h(j+1/2),    Stencil: [j+Ka, ......, j+Kb];    e.g.  Ka=-2, Kb=3 :  Stencil [j-2, j-1, j, j+1, j+2, j+3]
! Ka, Kb 是网格基架点的位置； 例如 Ka=-1, Kb=2 表示计算h(j+1/2)的基架点为 [j-1, j, j+1, j+2]
! v(0) ==>  v(j) ;  v(1) ==> v(j+1) ...
! hh ==> h(j+1/2)

    subroutine hh_weno5P(Ka, Kb, v, hh)          ! Ka=-2, Kb=2
       Use OCFD_constants
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
       real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
       real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

       S0 = 13.d0/12.d0*(v(0) - 2.d0*v(1) + v(2))**2 + 1.d0/4.d0*(3.d0*v(0) - 4.d0*v(1) + v(2))**2
       S1 = 13.d0/12.d0*(v(-1) - 2.d0*v(0) + v(1))**2 + 1.d0/4.d0*(v(-1) - v(1))**2
       S2 = 13.d0/12.d0*(v(-2) - 2.d0*v(-1) + v(0))**2 + 1.d0/4.d0*(v(-2) - 4.d0*v(-1) + 3.d0*v(0))**2
       a0 = C03/((ep + S0)**2)
       a1 = C13/((ep + S1)**2)
       a2 = C23/((ep + S2)**2)
       am = a0 + a1 + a2
       q03 = 1.d0/3.d0*v(0) + 5.d0/6.d0*v(1) - 1.d0/6.d0*v(2)
       q13 = -1.d0/6.d0*v(-1) + 5.d0/6.d0*v(0) + 1.d0/3.d0*v(1)
       q23 = 1.d0/3.d0*v(-2) - 7.d0/6.d0*v(-1) + 11.d0/6.d0*v(0)
       hh = (a0*q03 + a1*q13 + a2*q23)/am

    end

    subroutine hh_weno5M(Ka, Kb, v, hh)           ! Ka=-1,  Kb=3
       Use OCFD_constants
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
       real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
       real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0
       ! for flux-

       S0 = 13.d0/12.d0*(v(1) - 2.d0*v(0) + v(-1))**2 + 1.d0/4.d0*(3.d0*v(1) - 4.d0*v(0) + v(-1))**2
       S1 = 13.d0/12.d0*(v(2) - 2.d0*v(1) + v(0))**2 + 1.d0/4.d0*(v(2) - v(0))**2
       S2 = 13.d0/12.d0*(v(3) - 2.d0*v(2) + v(1))**2 + 1.d0/4.d0*(v(3) - 4.d0*v(2) + 3.d0*v(1))**2
       a0 = C03/((ep + S0)**2)
       a1 = C13/((ep + S1)**2)
       a2 = C23/((ep + S2)**2)
       am = a0 + a1 + a2
       q03 = 1.d0/3.d0*v(1) + 5.d0/6.d0*v(0) - 1.d0/6.d0*v(-1)
       q13 = -1.d0/6.d0*v(2) + 5.d0/6.d0*v(1) + 1.d0/3.d0*v(0)
       q23 = 1.d0/3.d0*v(3) - 7.d0/6.d0*v(2) + 11.d0/6.d0*v(1)
       hh = (a0*q03 + a1*q13 + a2*q23)/am

    end

!-----------------------------------------------------------------------------------------------------
!   7th order WENO-JS schemes

!----------------------------------------------------------------------------------
    subroutine hh_weno7P(Ka, Kb, v, hh)           ! Ka=-3,  Kb=3
       Use OCFD_constants
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
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

! 7th order WENO scheme
! 1  阶导数
       S10 = a11*v(-3) + a12*v(-2) + a13*v(-1) + a14*v(0)
       S11 = a21*v(-2) - v(-1) + a23*v(0) + a24*v(1)
       S12 = a31*v(-1) + a32*v(0) + v(1) + a34*v(2)
       S13 = a41*v(0) + a42*v(1) + a43*v(2) + a44*v(3)
       ! 2 阶导数
       S20 = -v(-3) + b12*v(-2) + b13*v(-1) + b14*v(0)
       S21 = v(-1) + b22*v(0) + v(1)
       S22 = v(0) + b22*v(1) + v(2)
       S23 = b41*v(0) + b42*v(1) + b43*v(2) - v(3)
! 3 阶导数
       S30 = -v(-3) + c12*(v(-2) - v(-1)) + v(0)
       S31 = -v(-2) + c12*(v(-1) - v(0)) + v(1)
       S32 = -v(-1) + c12*(v(0) - v(1)) + v(2)
       S33 = -v(0) + c12*(v(1) - v(2)) + v(3)

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
       q0 = e11*v(-3) + e12*v(-2) + e13*v(-1) + e14*v(0)
       q1 = e21*v(-2) + e22*v(-1) + e23*v(0) + e24*v(1)
       q2 = e31*v(-1) + e32*v(0) + e33*v(1) + e34*v(2)
       q3 = e41*v(0) + e42*v(1) + e43*v(2) + e44*v(3)

!  由4个4阶差分格式组合成1个7阶差分格式
!     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3
       hh = (a0*q0 + a1*q1 + a2*q2 + a3*q3)/am

    end

!----------WENO7 for flux-
    subroutine hh_weno7M(Ka, Kb, v, hh)          ! Ka=-2, Kb=4
       Use OCFD_constants
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
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

!      7th order WENO scheme
! 1  阶导数
       S10 = a11*v(4) + a12*v(3) + a13*v(2) + a14*v(1)
       S11 = a21*v(3) - v(2) + a23*v(1) + a24*v(0)
       S12 = a31*v(2) + a32*v(1) + v(0) + a34*v(-1)
       S13 = a41*v(1) + a42*v(0) + a43*v(-1) + a44*v(-2)
! 2 阶导数
       S20 = -v(4) + b12*v(3) + b13*v(2) + b14*v(1)
       S21 = v(2) + b22*v(1) + v(0)
       S22 = v(1) + b22*v(0) + v(-1)
       S23 = b41*v(1) + b42*v(0) + b43*v(-1) - v(-2)
! 3 阶导数
       S30 = -v(4) + c12*(v(3) - v(2)) + v(1)
       S31 = -v(3) + c12*(v(2) - v(1)) + v(0)
       S32 = -v(2) + c12*(v(1) - v(0)) + v(-1)
       S33 = -v(1) + c12*(v(0) - v(-1)) + v(-2)

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
       q0 = e11*v(4) + e12*v(3) + e13*v(2) + e14*v(1)
       q1 = e21*v(3) + e22*v(2) + e23*v(1) + e24*v(0)
       q2 = e31*v(2) + e32*v(1) + e33*v(0) + e34*v(-1)
       q3 = e41*v(1) + e42*v(0) + e43*v(-1) + e44*v(-2)

!  由4个4阶差分格式组合成1个7阶差分格式
       hh = (a0*q0 + a1*q1 + a2*q2 + a3*q3)/am

    end

!---------------------------------------------

!---------------------------------------------
! Optimized 6th order Monotonicity-Preserving Schemes
! Scheme by Leng Yan & Li Xinliang  see: Int. J. Numer. Meth. Fluids 2013; 73:560–577
!================================================================================
    subroutine hh_OMP6P(Ka, Kb, v, hh)      !Ka=-3, Kb=4
       Use OCFD_constants
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
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

       mid_nf = OMP_a1*v(4) + OMP_a2*v(3) + OMP_a3*v(2) + OMP_a4*v(1) + OMP_a5*v(0) + OMP_a6*v(-1) + OMP_a7*v(-2) + OMP_a8*v(-3)
       mp = v(0) + minmod2((v(1) - v(0)), kappa*(v(0) - v(-1)))
       if ((mid_nf - v(0))*(mid_nf - mp) .ge. ep) then
          d1 = v(-2) + v(0) - 2.d0*v(-1)
          d2 = v(-1) + v(1) - 2.d0*v(0)
          d3 = v(0) + v(2) - 2.d0*v(1)
          ul = v(0) + kappa*(v(0) - v(-1))
          md = 0.5d0*(v(0) + v(1)) - 0.5d0*minmod4(4.d0*d2 - d3, 4.d0*d3 - d2, d2, d3)
          lc = v(0) + 0.5d0*(v(0) - v(-1)) + kappa*minmod4(4.d0*d1 - d2, 4.d0*d2 - d1, d2, d1)/3.d0
          fmin = max(min(v(0), v(1), md), min(v(0), ul, lc))
          fmax = min(max(v(0), v(1), md), max(v(0), ul, lc))
          mid_nf = mid_nf + minmod2(fmax - mid_nf, fmin - mid_nf)
       end if
       hh = mid_nf

    end

    subroutine hh_OMP6M(Ka, Kb, v, hh)     ! Ka=-3, Kb=4
       Use OCFD_constants
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
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

       mid_nf = OMP_a1*v(-3) + OMP_a2*v(-2) + OMP_a3*v(-1) + OMP_a4*v(0) + OMP_a5*v(1) + OMP_a6*v(2) + OMP_a7*v(3) + OMP_a8*v(4)
       mp = v(1) + minmod2((v(0) - v(1)), kappa*(v(1) - v(2)))
       if ((mid_nf - v(1))*(mid_nf - mp) .ge. ep) then
          d1 = v(3) + v(1) - 2.d0*v(2)
          d2 = v(2) + v(0) - 2.d0*v(1)
          d3 = v(1) + v(-1) - 2.d0*v(0)
          !----
          ul = v(1) + kappa*(v(1) - v(2))
          md = 0.5d0*(v(1) + v(0)) - 0.5d0*minmod4(4.d0*d2 - d3, 4.d0*d3 - d2, d2, d3)
          lc = v(1) + 0.5d0*(v(1) - v(0)) + kappa*minmod4(4.d0*d1 - d2, 4.d0*d2 - d1, d2, d1)/3.d0
          fmin = max(min(v(1), v(0), md), min(v(1), ul, lc))
          fmax = min(max(v(1), v(0), md), max(v(1), ul, lc))
          mid_nf = mid_nf + minmod2(fmax - mid_nf, fmin - mid_nf)
       end if
       hh = mid_nf

    end

!------------------------------------------------------
! Boundary Scheme (WENO5 Delete Stencil type)
! 边界格式： 基于WENO5, 删除掉出界的模板

    subroutine hh_weno5P_boundary(Ka, Kb, v, hh, Del_Stencil)          ! Ka=-2, Kb=2
       Use OCFD_constants
       implicit none
       integer Ka, Kb, Del_Stencil
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
       real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
       real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0
!-------for flux+
       a0 = 0.d0; q03 = 0.d0
       a2 = 0.d0; q23 = 0.d0

       if (Del_Stencil .ne. DEL_RIGHT) then
          S0 = 13.d0/12.d0*(v(0) - 2.d0*v(1) + v(2))**2 + 1.d0/4.d0*(3.d0*v(0) - 4.d0*v(1) + v(2))**2
          a0 = C03/((ep + S0)**2)
          q03 = 1.d0/3.d0*v(0) + 5.d0/6.d0*v(1) - 1.d0/6.d0*v(2)
       end if

       S1 = 13.d0/12.d0*(v(-1) - 2.d0*v(0) + v(1))**2 + 1.d0/4.d0*(v(-1) - v(1))**2
       a1 = C13/((ep + S1)**2)
       q13 = -1.d0/6.d0*v(-1) + 5.d0/6.d0*v(0) + 1.d0/3.d0*v(1)

       if (Del_Stencil .ne. DEL_LIFT) then
          S2 = 13.d0/12.d0*(v(-2) - 2.d0*v(-1) + v(0))**2 + 1.d0/4.d0*(v(-2) - 4.d0*v(-1) + 3.d0*v(0))**2
          a2 = C23/((ep + S2)**2)
          q23 = 1.d0/3.d0*v(-2) - 7.d0/6.d0*v(-1) + 11.d0/6.d0*v(0)
       end if
       am = a0 + a1 + a2
       hh = (a0*q03 + a1*q13 + a2*q23)/am

    end

    subroutine hh_weno5M_boundary(Ka, Kb, v, hh, Del_Stencil)           ! Ka=-1,  Kb=3
       Use OCFD_constants
       implicit none
       integer Ka, Kb, Del_Stencil
       real(kind=OCFD_REAL_KIND)::  v(Ka:Kb), hh
       real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
       real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0
       ! for flux-
       a0 = 0.d0; q03 = 0.d0
       a2 = 0.d0; q23 = 0.d0

       if (Del_Stencil .ne. DEL_LIFT) then
          S0 = 13.d0/12.d0*(v(1) - 2.d0*v(0) + v(-1))**2 + 1.d0/4.d0*(3.d0*v(1) - 4.d0*v(0) + v(-1))**2
          a0 = C03/((ep + S0)**2)
          q03 = 1.d0/3.d0*v(1) + 5.d0/6.d0*v(0) - 1.d0/6.d0*v(-1)
       end if
       S1 = 13.d0/12.d0*(v(2) - 2.d0*v(1) + v(0))**2 + 1.d0/4.d0*(v(2) - v(0))**2
       a1 = C13/((ep + S1)**2)
       q13 = -1.d0/6.d0*v(2) + 5.d0/6.d0*v(1) + 1.d0/3.d0*v(0)

       if (Del_Stencil .ne. DEL_RIGHT) then
          S2 = 13.d0/12.d0*(v(3) - 2.d0*v(2) + v(1))**2 + 1.d0/4.d0*(v(3) - 4.d0*v(2) + 3.d0*v(1))**2
          a2 = C23/((ep + S2)**2)
          q23 = 1.d0/3.d0*v(3) - 7.d0/6.d0*v(2) + 11.d0/6.d0*v(1)
       end if
       am = a0 + a1 + a2
       hh = (a0*q03 + a1*q13 + a2*q23)/am

    end

!-------7th order low-dissipation scheme (mixing 7th upwind and 8th center)
    subroutine hh_UD7L_P(Ka, Kb, v, hh)
       use OCFD_constants
       use Scheme_para
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND)::  v(Ka:Kb), hh

       hh = Scheme%UD7L(1)*v(-3) + Scheme%UD7L(2)*v(-2) + Scheme%UD7L(3)*v(-1) + Scheme%UD7L(4)*v(0) &
            + Scheme%UD7L(5)*v(1) + Scheme%UD7L(6)*v(2) + Scheme%UD7L(7)*v(3) + Scheme%UD7L(8)*v(4)
    end

    subroutine hh_UD7L_M(Ka, Kb, v, hh)
       use OCFD_constants
       use Scheme_para
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND)::  v(Ka:Kb), hh

       hh = Scheme%UD7L(1)*v(4) + Scheme%UD7L(2)*v(3) + Scheme%UD7L(3)*v(2) + Scheme%UD7L(4)*v(1) &
            + Scheme%UD7L(5)*v(0) + Scheme%UD7L(6)*v(-1) + Scheme%UD7L(7)*v(-2) + Scheme%UD7L(8)*v(-3)
    end
!-----------------------------------------------------------------------------
! Hybrid scheme (UD7L + WENO7+ WENO5)
    subroutine hh_HybridP(Ka, Kb, v, hh, Ihybrid)          ! Ka=-4, Kb=4
       Use OCFD_constants
       use Scheme_para
       implicit none
       integer Ka, Kb, Ihybrid
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
       if (Ihybrid == 1) then       ! UD7L
          hh = Scheme%UD7L(1)*v(-3) + Scheme%UD7L(2)*v(-2) + Scheme%UD7L(3)*v(-1) + Scheme%UD7L(4)*v(0) &
               + Scheme%UD7L(5)*v(1) + Scheme%UD7L(6)*v(2) + Scheme%UD7L(7)*v(3) + Scheme%UD7L(8)*v(4)
       else if (Ihybrid == 2) then     ! WENO7 scheme
          call hh_weno7P(Ka, Kb, v, hh)
       else
          call hh_weno5P(Ka, Kb, v, hh)
       end if
    end

    subroutine hh_HybridM(Ka, Kb, v, hh, Ihybrid)           ! Ka=-4,  Kb=4
       Use OCFD_constants
       use Scheme_para
       implicit none
       integer Ka, Kb, Ihybrid
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
       real(kind=OCFD_REAL_KIND)::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23
       real(kind=OCFD_REAL_KIND), parameter::   ep = 1.d-6, C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0
       ! for flux-
       if (Ihybrid == 1) then
          hh = Scheme%UD7L(1)*v(4) + Scheme%UD7L(2)*v(3) + Scheme%UD7L(3)*v(2) + Scheme%UD7L(4)*v(1) &
               + Scheme%UD7L(5)*v(0) + Scheme%UD7L(6)*v(-1) + Scheme%UD7L(7)*v(-2) + Scheme%UD7L(8)*v(-3)
       else if (Ihybrid == 2) then      ! WENO7 scheme
          call hh_weno7M(Ka, Kb, v, hh)
       else
          call hh_weno5M(Ka, Kb, v, hh)
       end if
    end

!================================================================================

    subroutine hh_NND2P(Ka, Kb, v, hh)     ! Ka=-1, Kb=1
       Use OCFD_precision
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh, minmod2
       hh = v(0) + 0.5d0*minmod2(v(1) - v(0), v(0) - v(-1))     ! v(0) ==> v(j);  hh==> h(j+1/2)
    end

    subroutine hh_NND2M(Ka, Kb, v, hh)     ! Ka=0, Kb=2
       Use OCFD_precision
       implicit none
       integer Ka, Kb
       real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh, minmod2
       hh = v(1) - 0.5d0*minmod2(v(2) - v(1), v(1) - v(0))
    end
