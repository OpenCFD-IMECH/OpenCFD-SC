     subroutine Jac_Ghost_boundary
        use flow_data
        implicit none
        call Jac_Ghost_Extent_1st    ! 1st order exterpolation , in case Para%Ghost_Cell(1) .eq. 1
        call Jac_Ghost_Extent_2nd    ! 2nd order exterpolation , in case Para%Ghost_Cell(1) .eq. 1
        call Jac_Ghost_Userdef   ! User defined Ghost boundary , in case Para%Ghost_Cell(1) .eq. 1
     end

     subroutine Flow_Ghost_boundary
        use flow_data
        implicit none
        call Flow_Ghost_Extent_2nd  ! 2nd exterpolation, in case Scheme%Scheme_boundary(:) .eq. 1 or 2
        call Flow_Ghost_Userdef  ! User defined Ghost boundary, in case Scheme%Scheme_boundary(:) .eq. -2
     end

!--------------------------------------------------
! Ghost Cell for Jacobian coefficients ; 1st order exterpolation
     subroutine Jac_Ghost_Extent_1st
        use flow_data
        implicit none
        integer:: i, j, k
        if (npx .eq. 0 .and. Para%Ghost_Cell(1) .eq. 1) then  !i-
           do k = 1, nz
           do j = 1, ny
           do i = 1 - LAP, 0
              Akx(i, j, k) = Akx(1, j, k)
              Aky(i, j, k) = Aky(1, j, k)
              Akz(i, j, k) = Akz(1, j, k)
              Aix(i, j, k) = Aix(1, j, k)
              Aiy(i, j, k) = Aiy(1, j, k)
              Aiz(i, j, k) = Aiz(1, j, k)
              Asx(i, j, k) = Asx(1, j, k)
              Asy(i, j, k) = Asy(1, j, k)
              Asz(i, j, k) = Asz(1, j, k)
              Ajac(i, j, k) = Ajac(1, j, k)  ! A bug romoved, 2021-6-14

              Axx(i, j, k) = Axx(1, j, k) + (1.d0 - i)*(Axx(1, j, k) - Axx(2, j, k))
              Ayy(i, j, k) = Ayy(1, j, k) + (1.d0 - i)*(Ayy(1, j, k) - Ayy(2, j, k))
              Azz(i, j, k) = Azz(1, j, k) + (1.d0 - i)*(Azz(1, j, k) - Azz(2, j, k))
           end do
           end do
           end do
        end if

        if (npx .eq. npx0 - 1 .and. Para%Ghost_Cell(2) .eq. 1) then  ! i+
           do k = 1, nz
           do j = 1, ny
           do i = nx + 1, nx + LAP
              Akx(i, j, k) = Akx(nx, j, k)
              Aky(i, j, k) = Aky(nx, j, k)
              Akz(i, j, k) = Akz(nx, j, k)
              Aix(i, j, k) = Aix(nx, j, k)
              Aiy(i, j, k) = Aiy(nx, j, k)
              Aiz(i, j, k) = Aiz(nx, j, k)
              Asx(i, j, k) = Asx(nx, j, k)
              Asy(i, j, k) = Asy(nx, j, k)
              Asz(i, j, k) = Asz(nx, j, k)
              Ajac(i, j, k) = Ajac(nx, j, k)

              Axx(i, j, k) = Axx(nx, j, k) + (i - nx)*(Axx(nx, j, k) - Axx(nx - 1, j, k))
              Ayy(i, j, k) = Ayy(nx, j, k) + (i - nx)*(Ayy(nx, j, k) - Ayy(nx - 1, j, k))
              Azz(i, j, k) = Azz(nx, j, k) + (i - nx)*(Azz(nx, j, k) - Azz(nx - 1, j, k))
           end do
           end do
           end do

        end if

        if (npy .eq. 0 .and. Para%Ghost_Cell(3) .eq. 1) then      ! j-
           do k = 1, nz
           do j = 1 - LAP, 0
              do i = 1, nx
                 Akx(i, j, k) = Akx(i, 1, k)
                 Aky(i, j, k) = Aky(i, 1, k)
                 Akz(i, j, k) = Akz(i, 1, k)
                 Aix(i, j, k) = Aix(i, 1, k)
                 Aiy(i, j, k) = Aiy(i, 1, k)
                 Aiz(i, j, k) = Aiz(i, 1, k)
                 Asx(i, j, k) = Asx(i, 1, k)
                 Asy(i, j, k) = Asy(i, 1, k)
                 Asz(i, j, k) = Asz(i, 1, k)
                 Ajac(i, j, k) = Ajac(i, 1, k)

                 Axx(i, j, k) = Axx(i, 1, k) + (1.d0 - j)*(Axx(i, 1, k) - Axx(i, 2, k))
                 Ayy(i, j, k) = Ayy(i, 1, k) + (1.d0 - j)*(Ayy(i, 1, k) - Ayy(i, 2, k))
                 Azz(i, j, k) = Azz(i, 1, k) + (1.d0 - j)*(Azz(i, 1, k) - Azz(i, 2, k))
              end do
           end do
           end do
        end if

        if (npy .eq. npy0 - 1 .and. Para%Ghost_Cell(4) .eq. 1) then  ! j+
           do k = 1, nz
           do j = ny + 1, ny + LAP
           do i = 1, nx
              Akx(i, j, k) = Akx(i, ny, k)
              Aky(i, j, k) = Aky(i, ny, k)
              Akz(i, j, k) = Akz(i, ny, k)
              Aix(i, j, k) = Aix(i, ny, k)
              Aiy(i, j, k) = Aiy(i, ny, k)
              Aiz(i, j, k) = Aiz(i, ny, k)
              Asx(i, j, k) = Asx(i, ny, k)
              Asy(i, j, k) = Asy(i, ny, k)
              Asz(i, j, k) = Asz(i, ny, k)
              Ajac(i, j, k) = Ajac(i, ny, k)

              Axx(i, j, k) = Axx(i, ny, k) + (j - ny)*(Axx(i, ny, k) - Axx(i, ny - 1, k))
              Ayy(i, j, k) = Ayy(i, ny, k) + (j - ny)*(Ayy(i, ny, k) - Ayy(i, ny - 1, k))
              Azz(i, j, k) = Azz(i, ny, k) + (j - ny)*(Azz(i, ny, k) - Azz(i, ny - 1, k))
           end do
           end do
           end do
        end if

        if (npz .eq. 0 .and. Para%Ghost_Cell(5) .eq. 1) then      ! k-
           do k = 1 - LAP, 0
           do j = 1, ny
              do i = 1, nx
                 Akx(i, j, k) = Akx(i, j, 1)
                 Aky(i, j, k) = Aky(i, j, 1)
                 Akz(i, j, k) = Akz(i, j, 1)
                 Aix(i, j, k) = Aix(i, j, 1)
                 Aiy(i, j, k) = Aiy(i, j, 1)
                 Aiz(i, j, k) = Aiz(i, j, 1)
                 Asx(i, j, k) = Asx(i, j, 1)
                 Asy(i, j, k) = Asy(i, j, 1)
                 Asz(i, j, k) = Asz(i, j, 1)
                 Ajac(i, j, k) = Ajac(i, j, 1)

                 Axx(i, j, k) = Axx(i, j, 1) + (1.d0 - k)*(Axx(i, j, 1) - Axx(i, j, 2))
                 Ayy(i, j, k) = Ayy(i, j, 1) + (1.d0 - k)*(Ayy(i, j, 1) - Ayy(i, j, 2))
                 Azz(i, j, k) = Azz(i, j, 1) + (1.d0 - k)*(Azz(i, j, 1) - Azz(i, j, 2))
              end do
           end do
           end do
        end if

        if (npz .eq. npz0 - 1 .and. Para%Ghost_Cell(6) .eq. 1) then  ! k+
           do k = nz + 1, nz + LAP
           do j = 1, ny
           do i = 1, nx
              Akx(i, j, k) = Akx(i, j, nz)
              Aky(i, j, k) = Aky(i, j, nz)
              Akz(i, j, k) = Akz(i, j, nz)
              Aix(i, j, k) = Aix(i, j, nz)
              Aiy(i, j, k) = Aiy(i, j, nz)
              Aiz(i, j, k) = Aiz(i, j, nz)
              Asx(i, j, k) = Asx(i, j, nz)
              Asy(i, j, k) = Asy(i, j, nz)
              Asz(i, j, k) = Asz(i, j, nz)
              Ajac(i, j, k) = Ajac(i, j, nz)

              Axx(i, j, k) = Axx(i, j, nz) + (k - nz)*(Axx(i, j, nz) - Axx(i, j, nz - 1))
              Ayy(i, j, k) = Ayy(i, j, nz) + (k - nz)*(Ayy(i, j, nz) - Ayy(i, j, nz - 1))
              Azz(i, j, k) = Azz(i, j, nz) + (k - nz)*(Azz(i, j, nz) - Azz(i, j, nz - 1))
           end do
           end do
           end do
        end if
     end

!-----------------------------------------------
! 2nd Exter-ploation Ghost Cell , in case Para%Ghost_Cell(:) ==1 or 2
     subroutine Flow_Ghost_Extent_2nd
        use flow_data
        implicit none
        integer i, j, k
        Real(kind=OCFD_REAL_KIND), parameter:: s1 = 0.8d0, s2 = 1.2d0       ! limter

        if (npx == 0 .and. (Para%Ghost_Cell(1) == 1 .or. Para%Ghost_Cell(1) == 2)) then   ! i-
           do k = 1, nz
           do j = 1, ny
           do i = 1 - LAP, 0
              u(i, j, k) = 2.d0*u(1, j, k) - u(2 - i, j, k)
              v(i, j, k) = 2.d0*v(1, j, k) - v(2 - i, j, k)
              w(i, j, k) = 2.d0*w(1, j, k) - w(2 - i, j, k)

              d(i, j, k) = 2.d0*d(1, j, k) - d(2 - i, j, k)
              T(i, j, k) = 2.d0*T(1, j, k) - T(2 - i, j, k)
              if (d(i, j, k) < s1*d(1, j, k)) d(i, j, k) = s1*d(1, j, k)
              if (d(i, j, k) > s2*d(1, j, k)) d(i, j, k) = s2*d(1, j, k)
              if (T(i, j, k) < s1*T(1, j, k)) T(i, j, k) = s1*T(1, j, k)
              if (T(i, j, k) > s2*T(1, j, k)) T(i, j, k) = s2*T(1, j, k)

           end do
           end do
           end do
        end if

        if (npx == npx0 - 1 .and. (Para%Ghost_Cell(2) == 1 .or. Para%Ghost_Cell(2) == 2)) then   ! i+
           do k = 1, nz
           do j = 1, ny
           do i = nx + 1, nx + LAP
              u(i, j, k) = 2.d0*u(nx, j, k) - u(2*nx - i, j, k)
              v(i, j, k) = 2.d0*v(nx, j, k) - v(2*nx - i, j, k)
              w(i, j, k) = 2.d0*w(nx, j, k) - w(2*nx - i, j, k)

              d(i, j, k) = 2.d0*d(nx, j, k) - d(2*nx - i, j, k)
              T(i, j, k) = 2.d0*T(nx, j, k) - T(2*nx - i, j, k)

              if (d(i, j, k) < s1*d(nx, j, k)) d(i, j, k) = s1*d(nx, j, k)
              if (d(i, j, k) > s2*d(nx, j, k)) d(i, j, k) = s2*d(nx, j, k)
              if (T(i, j, k) < s1*T(nx, j, k)) T(i, j, k) = s1*T(nx, j, k)
              if (T(i, j, k) > s2*T(nx, j, k)) T(i, j, k) = s2*T(nx, j, k)

           end do
           end do
           end do
        end if

        if (npy == 0 .and. (Para%Ghost_Cell(3) == 1 .or. Para%Ghost_Cell(3) == 2)) then   ! j-
           do k = 1, nz
           do j = 1 - LAP, 0
           do i = 1, nx
              u(i, j, k) = 2.d0*u(i, 1, k) - u(i, 2 - j, k)
              v(i, j, k) = 2.d0*v(i, 1, k) - v(i, 2 - j, k)
              w(i, j, k) = 2.d0*w(i, 1, k) - w(i, 2 - j, k)

              d(i, j, k) = 2.d0*d(i, 1, k) - d(i, 2 - j, k)
              T(i, j, k) = 2.d0*T(i, 1, k) - T(i, 2 - j, k)
              if (d(i, j, k) < s1*d(i, 1, k)) d(i, j, k) = s1*d(i, 1, k)
              if (d(i, j, k) > s2*d(i, 1, k)) d(i, j, k) = s2*d(i, 1, k)
              if (T(i, j, k) < s1*T(i, 1, k)) T(i, j, k) = s1*T(i, 1, k)
              if (T(i, j, k) > s2*T(i, 1, k)) T(i, j, k) = s2*T(i, 1, k)
           end do
           end do
           end do
        end if

        if (npy == npy0 - 1 .and. (Para%Ghost_Cell(4) == 1 .or. Para%Ghost_Cell(4) == 2)) then   ! j+
           do k = 1, nz
           do j = ny + 1, ny + LAP
           do i = 1, nx
              u(i, j, k) = 2.d0*u(i, ny, k) - u(i, 2*ny - j, k)
              v(i, j, k) = 2.d0*v(i, ny, k) - v(i, 2*ny - j, k)
              w(i, j, k) = 2.d0*w(i, ny, k) - w(i, 2*ny - j, k)

              d(i, j, k) = 2.d0*d(i, ny, k) - d(i, 2*ny - j, k)
              T(i, j, k) = 2.d0*T(i, ny, k) - T(i, 2*ny - j, k)
              if (d(i, j, k) < s1*d(i, ny, k)) d(i, j, k) = s1*d(i, ny, k)
              if (d(i, j, k) > s2*d(i, ny, k)) d(i, j, k) = s2*d(i, ny, k)
              if (T(i, j, k) < s1*T(i, ny, k)) T(i, j, k) = s1*T(i, ny, k)
              if (T(i, j, k) > s2*T(i, ny, k)) T(i, j, k) = s2*T(i, ny, k)
           end do
           end do
           end do
        end if

        if (npz == 0 .and. (Para%Ghost_Cell(5) == 1 .or. Para%Ghost_Cell(5) == 2)) then   ! k-
           do k = 1 - LAP, 0
           do j = 1, ny
           do i = 1, nx
              u(i, j, k) = 2.d0*u(i, j, 1) - u(i, j, 2 - k)
              v(i, j, k) = 2.d0*v(i, j, 1) - v(i, j, 2 - k)
              w(i, j, k) = 2.d0*w(i, j, 1) - w(i, j, 2 - k)
              d(i, j, k) = 2.d0*d(i, j, 1) - d(i, j, 2 - k)
              T(i, j, k) = 2.d0*T(i, j, 1) - T(i, j, 2 - k)
              if (d(i, j, k) < s1*d(i, j, 1)) d(i, j, k) = s1*d(i, j, 1)
              if (d(i, j, k) > s2*d(i, j, 1)) d(i, j, k) = s2*d(i, j, 1)
              if (T(i, j, k) < s1*T(i, j, 1)) T(i, j, k) = s1*T(i, j, 1)
              if (T(i, j, k) > s2*T(i, j, 1)) T(i, j, k) = s2*T(i, j, 1)

           end do
           end do
           end do
        end if

        if (npz == npz0 - 1 .and. (Para%Ghost_Cell(6) == 1 .or. Para%Ghost_Cell(6) == 2)) then   ! k+
           do k = nz + 1, nz + LAP
           do j = 1, ny
           do i = 1, nx
              u(i, j, k) = 2.d0*u(i, j, nz) - u(i, j, 2*nz - k)
              v(i, j, k) = 2.d0*v(i, j, nz) - v(i, j, 2*nz - k)
              w(i, j, k) = 2.d0*w(i, j, nz) - w(i, j, 2*nz - k)
              d(i, j, k) = 2.d0*d(i, j, nz) - d(i, j, 2*nz - k)
              T(i, j, k) = 2.d0*T(i, j, nz) - T(i, j, 2*nz - k)
              if (d(i, j, k) < s1*d(i, j, nz)) d(i, j, k) = s1*d(i, j, nz)
              if (d(i, j, k) > s2*d(i, j, nz)) d(i, j, k) = s2*d(i, j, nz)
              if (T(i, j, k) < s1*T(i, j, nz)) T(i, j, k) = s1*T(i, j, nz)
              if (T(i, j, k) > s2*T(i, j, nz)) T(i, j, k) = s2*T(i, j, nz)

           end do
           end do
           end do
        end if
     end

!--------------------------------------------------
! Ghost Cell for Jacobian coefficients ; 2nd order exterpolation
     subroutine Jac_Ghost_Extent_2nd
        use flow_data
        implicit none
        integer:: i, j, k, i1, j1, k1
        if (npx .eq. 0 .and. Para%Ghost_Cell(1) .eq. 2) then  !i-
           do k = 1, nz
           do j = 1, ny
           do i = 1 - LAP, 0
              i1 = 2 - i
              Axx(i, j, k) = 2.d0*Axx(1, j, k) - Axx(i1, j, k)
              Ayy(i, j, k) = 2.d0*Ayy(1, j, k) - Ayy(i1, j, k)
              Azz(i, j, k) = 2.d0*Azz(1, j, k) - Azz(i1, j, k)
           end do
           end do
           end do

           call comput_Jacobian3d_Ghost(1)           ! Comput Jocabian coefficient of the Ghost Cells

        end if

        if (npx .eq. npx0 - 1 .and. Para%Ghost_Cell(2) .eq. 2) then  ! i+
           do k = 1, nz
           do j = 1, ny
           do i = nx + 1, nx + LAP
              i1 = 2*nx - i
              Axx(i, j, k) = 2.d0*Axx(nx, j, k) - Axx(i1, j, k)
              Ayy(i, j, k) = 2.d0*Ayy(nx, j, k) - Ayy(i1, j, k)
              Azz(i, j, k) = 2.d0*Azz(nx, j, k) - Azz(i1, j, k)
           end do
           end do
           end do

           call comput_Jacobian3d_Ghost(2)           ! Comput Jocabian coefficient of the Ghost Cells
        end if

        if (npy .eq. 0 .and. Para%Ghost_Cell(3) .eq. 2) then      ! j-
           do k = 1, nz
           do j = 1 - LAP, 0
              do i = 1, nx
                 j1 = 2 - j
                 Axx(i, j, k) = 2.d0*Axx(i, 1, k) - Axx(i, j1, k)
                 Ayy(i, j, k) = 2.d0*Ayy(i, 1, k) - Ayy(i, j1, k)
                 Azz(i, j, k) = 2.d0*Azz(i, 1, k) - Azz(i, j1, k)
              end do
           end do
           end do
           call comput_Jacobian3d_Ghost(3)           ! Comput Jocabian coefficient of the Ghost Cells
        end if

        if (npy .eq. npy0 - 1 .and. Para%Ghost_Cell(4) .eq. 2) then  ! j+
           do k = 1, nz
           do j = ny + 1, ny + LAP
           do i = 1, nx
              j1 = 2*ny - j
              Axx(i, j, k) = 2.d0*Axx(i, ny, k) - Axx(i, j1, k)
              Ayy(i, j, k) = 2.d0*Ayy(i, ny, k) - Ayy(i, j1, k)
              Azz(i, j, k) = 2.d0*Azz(i, ny, k) - Azz(i, j1, k)
           end do
           end do
           end do
           call comput_Jacobian3d_Ghost(4)           ! Comput Jocabian coefficient of the Ghost Cells
        end if

        if (npz .eq. 0 .and. Para%Ghost_Cell(5) .eq. 2) then      ! k-
           do k = 1 - LAP, 0
           do j = 1, ny
              do i = 1, nx
                 k1 = 2 - k
                 Axx(i, j, k) = 2.d0*Axx(i, j, 1) - Axx(i, j, k1)
                 Ayy(i, j, k) = 2.d0*Ayy(i, j, 1) - Ayy(i, j, k1)
                 Azz(i, j, k) = 2.d0*Azz(i, j, 1) - Azz(i, j, k1)
              end do
           end do
           end do
           call comput_Jacobian3d_Ghost(5)           ! Comput Jocabian coefficient of the Ghost Cells
        end if

        if (npz .eq. npz0 - 1 .and. Para%Ghost_Cell(6) .eq. 2) then  ! k+
           do k = nz + 1, nz + LAP
           do j = 1, ny
           do i = 1, nx
              k1 = 2*nz - k
              Axx(i, j, k) = 2.d0*Axx(i, j, nz) - Axx(i, j, k1)
              Ayy(i, j, k) = 2.d0*Ayy(i, j, nz) - Ayy(i, j, k1)
              Azz(i, j, k) = 2.d0*Azz(i, j, nz) - Azz(i, j, k1)
           end do
           end do
           end do
           call comput_Jacobian3d_Ghost(6)           ! Comput Jocabian coefficient of the Ghost Cells
        end if

!    test test ------------------------
!        if(npx .eq. 11 .and. npy .eq. 0 .and. npz .eq. 4) then
!                open(99,file="test-jac.dat")
!                do j=1-LAP, 5
!                write(99,"(15F16.8)") 1.d0*j, Axx(1,j,1),Ayy(1,j,1),Azz(1,j,1),Akx(1,j,1),Aky(1,j,1),Akz(1,j,1), &
!                      Aix(1,j,1),Aiy(1,j,1),Aiz(1,j,1),Asx(1,j,1),Asy(1,j,1),Asz(1,j,1),Ajac(1,j,1)
!            enddo
!                close(99)
!                endif

     end

!--------Comput Jocabian coefficients at Ghost Cells
     subroutine comput_Jacobian3d_Ghost(nb)
        Use flow_data
        implicit none
        integer:: nb, ib, ie, jb, je, kb, ke         ! Jocabian data range
        integer::  ib1, ie1, jb1, je1, kb1, ke1               ! coordinate data range
        real(kind=OCFD_REAL_KIND):: xi1, xj1, xk1, yi1, yj1, yk1, zi1, zj1, zk1, Jac1
        integer:: i, j, k

        if (nb .eq. 1) then            ! i-
           ib = 1 - LAP; ie = 0; jb = 1; je = ny; kb = 1; ke = nz
           ib1 = ib; ie1 = 1; jb1 = jb; je1 = je; kb1 = kb; ke1 = ke
        else if (nb .eq. 2) then      ! i+
           ib = nx + 1; ie = nx + LAP; jb = 1; je = ny; kb = 1; ke = nz
           ib1 = nx; ie1 = ie; jb1 = jb; je1 = je; kb1 = kb; ke1 = ke
        else if (nb .eq. 3) then
           ib = 1; ie = nx; jb = 1 - LAP; je = 0; kb = 1; ke = nz
           ib1 = ib; ie1 = ie; jb1 = jb; je1 = 1; kb1 = kb; ke1 = ke
        else if (nb .eq. 4) then
           ib = 1; ie = nx; jb = ny + 1; je = ny + LAP; kb = 1; ke = nz
           ib1 = ib; ie1 = ie; jb1 = ny; je1 = je; kb1 = kb; ke1 = ke
        else if (nb .eq. 5) then
           ib = 1; ie = nx; jb = 1; je = ny; kb = 1 - LAP; ke = 0
           ib1 = ib; ie1 = ie; jb1 = jb; je1 = je; kb1 = kb; ke1 = 1
        else if (nb .eq. 6) then
           ib = 1; ie = nx; jb = 1; je = ny; kb = nz + 1; ke = nz + LAP
           ib1 = ib; ie1 = ie; jb1 = jb; je1 = je; kb1 = nz; ke1 = ke
        end if

        do k = kb, ke
        do j = jb, je
        do i = ib, ie

           if (i .eq. ib1) then
              xi1 = (-3.d0*Axx(i, j, k) + 4.d0*Axx(i + 1, j, k) - Axx(i + 2, j, k))/(2.d0*hx)     ! 2nd one-side scheme
              yi1 = (-3.d0*Ayy(i, j, k) + 4.d0*Ayy(i + 1, j, k) - Ayy(i + 2, j, k))/(2.d0*hx)
              zi1 = (-3.d0*Azz(i, j, k) + 4.d0*Azz(i + 1, j, k) - Azz(i + 2, j, k))/(2.d0*hx)
           else if (i .eq. ie1) then
              xi1 = (Axx(i - 2, j, k) - 4.d0*Axx(i - 1, j, k) + 3.d0*Axx(i, j, k))/(2.d0*hx)  ! 2nd one-side scheme
              yi1 = (Ayy(i - 2, j, k) - 4.d0*Ayy(i - 1, j, k) + 3.d0*Ayy(i, j, k))/(2.d0*hx)  ! 2nd one-side scheme
              zi1 = (Azz(i - 2, j, k) - 4.d0*Azz(i - 1, j, k) + 3.d0*Azz(i, j, k))/(2.d0*hx)  ! 2nd one-side scheme
           else if (i .eq. ib1 + 1 .or. i .eq. ie1 - 1) then
              xi1 = (Axx(i + 1, j, k) - Axx(i - 1, j, k))/(2.d0*hx)                             ! 2nd centeral scheme
              yi1 = (Ayy(i + 1, j, k) - Ayy(i - 1, j, k))/(2.d0*hx)
              zi1 = (Azz(i + 1, j, k) - Azz(i - 1, j, k))/(2.d0*hx)
           else
              xi1 = (8.d0*(Axx(i + 1, j, k) - Axx(i - 1, j, k)) - (Axx(i + 2, j, k) - Axx(i - 2, j, k)))/(12.d0*hx)   ! 4th central scheme
              yi1 = (8.d0*(Ayy(i + 1, j, k) - Ayy(i - 1, j, k)) - (Ayy(i + 2, j, k) - Ayy(i - 2, j, k)))/(12.d0*hx)
              zi1 = (8.d0*(Azz(i + 1, j, k) - Azz(i - 1, j, k)) - (Azz(i + 2, j, k) - Azz(i - 2, j, k)))/(12.d0*hx)
           end if

           if (j .eq. jb1) then
              xj1 = (-3.d0*Axx(i, j, k) + 4.d0*Axx(i, j + 1, k) - Axx(i, j + 2, k))/(2.d0*hy)     ! 2nd one-side scheme
              yj1 = (-3.d0*Ayy(i, j, k) + 4.d0*Ayy(i, j + 1, k) - Ayy(i, j + 2, k))/(2.d0*hy)
              zj1 = (-3.d0*Azz(i, j, k) + 4.d0*Azz(i, j + 1, k) - Azz(i, j + 2, k))/(2.d0*hy)

           else if (j .eq. je1) then
              xj1 = (Axx(i, j - 2, k) - 4.d0*Axx(i, j - 1, k) + 3.d0*Axx(i, j, k))/(2.d0*hy)  ! 2nd one-side scheme
              yj1 = (Ayy(i, j - 2, k) - 4.d0*Ayy(i, j - 1, k) + 3.d0*Ayy(i, j, k))/(2.d0*hy)
              zj1 = (Azz(i, j - 2, k) - 4.d0*Azz(i, j - 1, k) + 3.d0*Azz(i, j, k))/(2.d0*hy)
           else if (j .eq. jb1 + 1 .or. j .eq. je1 - 1) then
              xj1 = (Axx(i, j + 1, k) - Axx(i, j - 1, k))/(2.d0*hy)                      ! 2nd centeral scheme
              yj1 = (Ayy(i, j + 1, k) - Ayy(i, j - 1, k))/(2.d0*hy)
              zj1 = (Azz(i, j + 1, k) - Azz(i, j - 1, k))/(2.d0*hy)
           else
              xj1 = (8.d0*(Axx(i, j + 1, k) - Axx(i, j - 1, k)) - (Axx(i, j + 2, k) - Axx(i, j - 2, k)))/(12.d0*hy)   ! 4th central scheme
              yj1 = (8.d0*(Ayy(i, j + 1, k) - Ayy(i, j - 1, k)) - (Ayy(i, j + 2, k) - Ayy(i, j - 2, k)))/(12.d0*hy)
              zj1 = (8.d0*(Azz(i, j + 1, k) - Azz(i, j - 1, k)) - (Azz(i, j + 2, k) - Azz(i, j - 2, k)))/(12.d0*hy)
           end if

           if (k .eq. kb1) then
              xk1 = (-3.d0*Axx(i, j, k) + 4.d0*Axx(i, j, k + 1) - Axx(i, j, k + 2))/(2.d0*hz)     ! 2nd one-side scheme
              yk1 = (-3.d0*Ayy(i, j, k) + 4.d0*Ayy(i, j, k + 1) - Ayy(i, j, k + 2))/(2.d0*hz)
              zk1 = (-3.d0*Azz(i, j, k) + 4.d0*Azz(i, j, k + 1) - Azz(i, j, k + 2))/(2.d0*hz)
           else if (k .eq. ke1) then
              xk1 = (Axx(i, j, k - 2) - 4.d0*Axx(i, j, k - 1) + 3.d0*Axx(i, j, k))/(2.d0*hz)  ! 2nd one-side scheme
              yk1 = (Ayy(i, j, k - 2) - 4.d0*Ayy(i, j, k - 1) + 3.d0*Ayy(i, j, k))/(2.d0*hz)
              zk1 = (Azz(i, j, k - 2) - 4.d0*Azz(i, j, k - 1) + 3.d0*Azz(i, j, k))/(2.d0*hz)
           else if (k .eq. kb1 + 1 .or. k .eq. ke1 - 1) then
              xk1 = (Axx(i, j, k + 1) - Axx(i, j, k - 1))/(2.d0*hz)                             ! 2nd centeral scheme
              yk1 = (Ayy(i, j, k + 1) - Ayy(i, j, k - 1))/(2.d0*hz)
              zk1 = (Azz(i, j, k + 1) - Azz(i, j, k - 1))/(2.d0*hz)
           else
              xk1 = (8.d0*(Axx(i, j, k + 1) - Axx(i, j, k - 1)) - (Axx(i, j, k + 2) - Axx(i, j, k - 2)))/(12.d0*hz)   ! 4th central scheme
              yk1 = (8.d0*(Ayy(i, j, k + 1) - Ayy(i, j, k - 1)) - (Ayy(i, j, k + 2) - Ayy(i, j, k - 2)))/(12.d0*hz)
              zk1 = (8.d0*(Azz(i, j, k + 1) - Azz(i, j, k - 1)) - (Azz(i, j, k + 2) - Azz(i, j, k - 2)))/(12.d0*hz)
           end if

           Jac1 = 1.d0/(xi1*yj1*zk1 + yi1*zj1*xk1 + zi1*xj1*yk1 - zi1*yj1*xk1 - yi1*xj1*zk1 - xi1*zj1*yk1)   ! 1./Jocabian = d(x,y,z)/d(i,j,k)
           Ajac(i, j, k) = Jac1
           Akx(i, j, k) = Jac1*(yj1*zk1 - zj1*yk1)
           Aky(i, j, k) = Jac1*(zj1*xk1 - xj1*zk1)
           Akz(i, j, k) = Jac1*(xj1*yk1 - yj1*xk1)
           Aix(i, j, k) = Jac1*(yk1*zi1 - zk1*yi1)
           Aiy(i, j, k) = Jac1*(zk1*xi1 - xk1*zi1)
           Aiz(i, j, k) = Jac1*(xk1*yi1 - yk1*xi1)
           Asx(i, j, k) = Jac1*(yi1*zj1 - zi1*yj1)
           Asy(i, j, k) = Jac1*(zi1*xj1 - xi1*zj1)
           Asz(i, j, k) = Jac1*(xi1*yj1 - yi1*xj1)
           if (Jac1 .lt. 0) then
              print *, " Jocabian < 0 !!! , Jac=", Jac1
              print *, "i,j,k=", i_offset(npx) + i - 1, j_offset(npy) + j - 1, k_offset(npz) + k - 1
           end if
        end do
        end do
        end do

     end

