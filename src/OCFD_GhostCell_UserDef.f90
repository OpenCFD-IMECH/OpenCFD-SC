! User Defined Ghost-Cell boundary, for Jacobian and flow
!--------------------------------------------------
! Ghost Cell for Jacobian coefficients (User defined)
     subroutine Jac_Ghost_Userdef
        use flow_data
        implicit none
        if (Para%IBC .eq. BC_SweptCorner) then
           call Jac_Ghost_Sweptcorner
        else if (Para%IBC .eq. BC_BoundaryLayer) then
           call Jac_Ghost_flatplate                    ! only work in flatplate case  ( invalid  for compression corner or curve wall case)
        end if
     end

     subroutine Flow_Ghost_Userdef
        use flow_data
        implicit none
        if (Para%IBC .eq. BC_SweptCorner) then
           call Flow_Ghost_Sweptcorner
        else if (Para%IBC .eq. BC_BoundaryLayer) then
           call Flow_Ghost_flatplate
        end if
     end

! Ghost Cell for Jacobian coefficients for sweptcorner (symmetry in plane z=0)
     subroutine Jac_Ghost_Sweptcorner
        use flow_data
        implicit none
        integer:: i, j, k, k1
        if (npz .eq. 0 .and. Para%Ghost_Cell(5) .eq. -1) then  !k-
           do k = 1 - LAP, 0
              k1 = 1 - k
              do j = 1, ny
              do i = 1, nx
                 Akx(i, j, k) = Akx(i, j, k1)
                 Aky(i, j, k) = Aky(i, j, k1)
                 Akz(i, j, k) = -Akz(i, j, k1)
                 Aix(i, j, k) = Aix(i, j, k1)
                 Aiy(i, j, k) = Aiy(i, j, k1)
                 Aiz(i, j, k) = -Aiz(i, j, k1)
                 Asx(i, j, k) = -Asx(i, j, k1)
                 Asy(i, j, k) = -Asy(i, j, k1)
                 Asz(i, j, k) = Asz(i, j, k1)
                 Ajac(i, j, k) = Ajac(i, j, k1)

                 Axx(i, j, k) = Axx(i, j, k1)
                 Ayy(i, j, k) = Ayy(i, j, k1)
                 Azz(i, j, k) = -Azz(i, j, k1)
              end do
              end do
           end do
        end if
     end

!-----------------------------------------------
!  Ghost Cell for Sweptcorner (z=0 : symmetry )
     subroutine Flow_Ghost_Sweptcorner
        use flow_data
        implicit none
        integer i, j, k, k1

        if (npz == 0 .and. Para%Ghost_Cell(5) .eq. -1) then   ! k-
           do k = 1 - LAP, 0
              k1 = 1 - k
              do j = 1, ny
              do i = 1, nx
                 u(i, j, k) = u(i, j, k1)
                 v(i, j, k) = v(i, j, k1)
                 w(i, j, k) = -w(i, j, k1)
                 d(i, j, k) = d(i, j, k1)
                 T(i, j, k) = T(i, j, k1)
              end do
              end do
           end do
        end if

     end

! Used for flatplate boundary layers at y=0  (symmetry in y=0)
     subroutine Jac_Ghost_flatplate
        use flow_data
        implicit none
        integer:: i, j, k, j1
        if (npy .eq. 0 .and. Para%Ghost_Cell(3) .eq. -1) then  !j-
           do k = 1, nz
           do j = 1 - LAP, 0
              j1 = 2 - j                    !  0 <--->  2 ;   -1 <----> 3
              do i = 1, nx
                 Akx(i, j, k) = Akx(i, j1, k)
                 Aky(i, j, k) = -Aky(i, j1, k)
                 Akz(i, j, k) = Akz(i, j1, k)
                 Aix(i, j, k) = -Aix(i, j1, k)
                 Aiy(i, j, k) = Aiy(i, j1, k)
                 Aiz(i, j, k) = -Aiz(i, j1, k)
                 Asx(i, j, k) = Asx(i, j1, k)
                 Asy(i, j, k) = -Asy(i, j1, k)
                 Asz(i, j, k) = Asz(i, j1, k)
                 Ajac(i, j, k) = Ajac(i, j1, k)

                 Axx(i, j, k) = Axx(i, j1, k)
                 Ayy(i, j, k) = -Ayy(i, j1, k)
                 Azz(i, j, k) = Azz(i, j1, k)
              end do
           end do
           end do
        end if
     end

!-----------------------------------------------
!  Ghost Cell for flatplate (y=0 : symmetry )
     subroutine Flow_Ghost_flatplate
        use flow_data
        implicit none
        integer i, j, k, j1
        Real(kind=OCFD_REAL_KIND), parameter:: s1 = 0.5d0, s2 = 2.d0       ! limter

        if (npy == 0 .and. Para%Ghost_Cell(3) == -1) then   ! j-

           do k = 1, nz
           do j = 1 - LAP, 0
              j1 = 2 - j
              do i = 1, nx
                 u(i, j, k) = 2.d0*u(i, 1, k) - u(i, j1, k)
                 v(i, j, k) = 2.d0*v(i, 1, k) - v(i, j1, k)
                 w(i, j, k) = 2.d0*w(i, 1, k) - w(i, j1, k)

                 d(i, j, k) = 2.d0*d(i, 1, k) - d(i, j1, k)
                 T(i, j, k) = 2.d0*T(i, 1, k) - T(i, j1, k)
                 if (d(i, j, k) < s1*d(i, 1, k)) d(i, j, k) = s1*d(i, 1, k)
                 if (d(i, j, k) > s2*d(i, 1, k)) d(i, j, k) = s2*d(i, 1, k)
                 if (T(i, j, k) < s1*T(i, 1, k)) T(i, j, k) = s1*T(i, 1, k)
                 if (T(i, j, k) > s2*T(i, 1, k)) T(i, j, k) = s2*T(i, 1, k)
              end do
           end do
           end do
        end if
     end
