! -------------------------------------------------------------------------------
! read flow parameters  (Namelist type)
         subroutine read_parameter
            use flow_para
            implicit none
            integer::  IF_Scheme_Character, Iperiodic_X, Iperiodic_Y, Iperiodic_Z, IF_Viscous, Istep_Show, Istep_Save, &
                      Iflag_Gridtype, Scheme_boundary(6), IF_Mass_Force, ANA_Number, NFiltering, Flux_Splitting, &
                      Ghost_Cell(6), ierr
            character(len=50):: Scheme_Invis, Scheme_Vis, BoundaryCondition
            real(kind=OCFD_REAL_KIND):: Ma, Re, gamma, Pr, Ref_Amu_T0, dt, End_time, AoA, &
                                        Periodic_ISpan(3), Periodic_JSpan(3), Periodic_KSpan(3), &
                                        BC_Para(100), Mass_Force(3), ANA_Para(100, 10), Filter_Para(100, 10), &
                                        Hybrid_para(100), UD7L_Diss
            namelist /control_opencfd/ Ma, Re, gamma, Pr, Ref_Amu_T0, dt, End_time, AoA, &
               nx_global, ny_global, nz_global, npx0, npy0, npz0, LAP, &
               Iperiodic_X, Iperiodic_Y, Iperiodic_Z, &
               Scheme_Invis, Scheme_Vis, IF_Scheme_Character, IF_Viscous, &
               Istep_Show, Istep_Save, Iflag_Gridtype, &
               Periodic_ISpan, Periodic_JSpan, Periodic_KSpan, &
               BoundaryCondition, BC_Para, Scheme_boundary, IF_Mass_Force, Mass_Force, &
               ANA_Number, ANA_Para, NFiltering, Filter_Para, UD7L_Diss, Hybrid_para, Flux_Splitting, &
               Ghost_Cell

            integer:: nparameters(100)
            real(kind=OCFD_REAL_KIND):: rparameters(100)
            real*8, parameter::  PI = 3.1415926535897932d0
!--------Set defalut parameter ----------------
            gamma = 1.4d0
            Re = 100.d0
            Ma = 1.d0
            Pr = 0.7d0
            npx0 = 1
            npy0 = 1
            npz0 = 1
            LAP = 4
            Iperiodic_X = 0
            Iperiodic_Y = 0
            Iperiodic_Z = 0
            Scheme_Invis = "WENO5"
            Scheme_Vis = "CD6"
            IF_Scheme_Character = 0
            Ref_Amu_T0 = 288.15d0
            dt = 1.d-4
            End_time = 100.d0
            Iflag_Gridtype = GRID3D
            Periodic_ISpan(:) = 0.d0
            Periodic_JSpan(:) = 0.d0
            Periodic_KSpan(:) = 0.d0
            BC_Para(:) = 0.d0
            BoundaryCondition = "BC_None"
            Istep_show = 1
            Istep_save = 1000
            IF_Viscous = 1
            Scheme_boundary(:) = 0
            IF_Mass_Force = 0
            Mass_Force(:) = 0.d0
            AoA = 0.d0
            ANA_Number = 0
            ANA_Para = 0.d0
            NFiltering = 0
            Filter_Para = 0.d0
            UD7L_Diss = 1.d0    ! Disspation for Low-dissipative Upwind difference scheme (0-1,   0 CD8 scheme, 1 UD7 scheme)
            Hybrid_para(:) = 1.d0  ! parameters for Hybrid scheme
            Flux_Splitting = OCFD_Split_SW
            Ghost_Cell(:) = 0

!----------------------------------------------

            if (my_id .eq. 0) then
               print *, "---------OpenCFD-SC Version 2-------------"
               print *, "-       code by Li Xinliang              -"
               print *, "-    LHD, Institue of Mechanics, CAS     -"
               print *, "-      lixl@imech.ac.cn                  -"
               print *, "------------------------------------------"

               open (99, file="opencfd2.in")
               read (99, nml=control_opencfd)
               close (99)
               print *, "Re=", Re
               print *, "Ma=", Ma
               print *, "Scheme_Invis=", Scheme_Invis
               print *, "Scheme_Vis=", Scheme_Vis
               print *, "Flux_Splitting (1 SW; 2 LLF)=  ", Flux_Splitting
               print *, "IF_Scheme_Character (0/1) = ", IF_Scheme_Character
               print *, "-----------------------------------------------------"

               call capitalize(Scheme_Invis)
               call capitalize(Scheme_Vis)
               call capitalize(BoundaryCondition)

               call find_scheme(Scheme_Invis, Scheme%Scheme_Invis)     ! Scheme for inviscous terms
               call find_scheme(Scheme_Vis, Scheme%Scheme_Vis)             ! Scheme for viscous terms
               call find_bc(BoundaryCondition, Para%IBC)           ! boundary condition
               if (ANA_Number > 10 .or. NFiltering > 10) then
                  print *, "ANA_Number and NFiltering should <=10, please check !!!"
                  stop
               end if

               Para%BC_Para(1:100) = BC_Para(1:100)
               Para%ANA_Para(:, :) = ANA_Para(:, :)
               Para%Filter(:, :) = Filter_Para(:, :)
               Scheme%Hybrid_para(:) = Hybrid_para(:)

               nparameters(1) = nx_global
               nparameters(2) = ny_global
               nparameters(3) = nz_global

               nparameters(4) = npx0
               nparameters(5) = npy0
               nparameters(6) = npz0

               nparameters(7) = LAP
               nparameters(8) = Iperiodic_X
               nparameters(9) = Iperiodic_Y
               nparameters(10) = Iperiodic_Z

               nparameters(11) = Scheme%Scheme_Invis
               nparameters(12) = Scheme%Scheme_Vis
               nparameters(13) = If_Scheme_Character
               nparameters(14) = IF_Viscous
               nparameters(15) = Istep_Show
               nparameters(16) = Istep_Save
               nparameters(17) = Iflag_Gridtype
               nparameters(18) = Para%IBC
               nparameters(19:24) = Scheme_boundary(1:6)
               nparameters(25) = IF_Mass_Force
               nparameters(26) = ANA_Number
               nparameters(27) = NFiltering
               nparameters(28) = Flux_Splitting
               nparameters(29:34) = Ghost_Cell(1:6)

               rparameters(1) = Ma
               rparameters(2) = Re
               rparameters(3) = gamma
               rparameters(4) = Pr
               rparameters(5) = Ref_Amu_T0
               rparameters(6) = dt
               rparameters(7) = End_time
               rparameters(8:10) = Periodic_ISpan
               rparameters(11:13) = Periodic_JSpan
               rparameters(14:16) = Periodic_KSpan
               rparameters(17:19) = Mass_Force
               rparameters(20) = AoA*PI/180.d0   ! angle of attack (in degree)
               rparameters(21) = UD7L_Diss
            end if

! Broadcast the parameters to all MPI procs

            call MPI_bcast(nparameters(1), 100, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(rparameters(1), 100, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)

            nx_global = nparameters(1)
            ny_global = nparameters(2)
            nz_global = nparameters(3)

            npx0 = nparameters(4)
            npy0 = nparameters(5)
            npz0 = nparameters(6)

            LAP = nparameters(7)
            Para%Iperiodic_X = nparameters(8)
            Para%Iperiodic_Y = nparameters(9)
            Para%Iperiodic_Z = nparameters(10)

            Scheme%Scheme_Invis = nparameters(11)
            Scheme%Scheme_Vis = nparameters(12)
            Para%IF_Scheme_Character = nparameters(13)
            Para%IF_Viscous = nparameters(14)
            Para%Istep_Show = nparameters(15)
            Para%Istep_Save = nparameters(16)
            Para%Iflag_Gridtype = nparameters(17)
            Para%IBC = nparameters(18)
            Scheme%Scheme_boundary(1:6) = nparameters(19:24)
            Para%IF_Mass_Force = nparameters(25)
            Para%ANA_Number = nparameters(26)
            Para%NFiltering = nparameters(27)
            Para%Flux_Splitting = nparameters(28)
            Para%Ghost_Cell(1:6) = nparameters(29:34)

!----------------------------

            Para%Ma = rparameters(1)
            Para%Re = rparameters(2)
            Para%gamma = rparameters(3)
            Para%Pr = rparameters(4)
            Para%Ref_Amu_T0 = rparameters(5)
            Para%dt = rparameters(6)
            Para%End_time = rparameters(7)
            Para%Periodic_ISpan(:) = rparameters(8:10)
            Para%Periodic_JSpan(:) = rparameters(11:13)
            Para%Periodic_KSpan(:) = rparameters(14:16)
            Para%Mass_Force(:) = rparameters(17:19)
            Para%AoA = rparameters(20)
            Scheme%UD7L_Diss = rparameters(21)

            call MPI_bcast(Para%BC_Para(1), 100, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(Para%ANA_Para(1, 1), 1000, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(Para%Filter(1, 1), 1000, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(Scheme%Hybrid_para(1), 100, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)
         end

!-----------------------------------------------
         subroutine find_scheme(Schm1, scheme)
            use OCFD_constants
            implicit none
            character(len=50):: Schm1
            integer:: scheme

            select case (trim(Schm1))

            case ("WENO5")
               Scheme = OCFD_Scheme_WENO5
            case ("WENO7")
               Scheme = OCFD_Scheme_WENO7
            case ("OMP6")
               Scheme = OCFD_Scheme_OMP6
            case ("CD6")
               Scheme = OCFD_Scheme_CD6
            case ("CD8")
               Scheme = OCFD_Scheme_CD8
            case ("UD7L")
               Scheme = OCFD_Scheme_UD7L
            case ("HYBRID")
               scheme = OCFD_Scheme_Hybrid
            case ("SCHEME_USER")
               Scheme = OCFD_Scheme_USER
            case default
               print *, "This Scheme  is not supported !!", trim(Schm1)
               print *, "Scheme_Invis can be WENO5, WENO7 or OMP6 "
               print *, "Scheme_Vis  can be CD6 or CD8  (6th or 8th Centeral Scheme)"
               stop
            end select
         end

         subroutine find_bc(bc, ibc)
            use OCFD_constants
            implicit none
            character(len=50):: bc
            integer:: ibc

            select case (trim(bc))

            case ("BC_NONE")
               ibc = BC_None
            case ("BC_BLUNT2D")
               ibc = BC_Blunt2d
            case ("BC_BOUNDARYLAYER")
               ibc = BC_BoundaryLayer
            case ("BC_SWEPTCORNER")
               ibc = BC_SweptCorner
            case ("BC_USER_DEF")
               ibc = BC_User_Def
            case default
               print *, "This Boundary Condition  is not supported !!", trim(bc)
               print *, "Please read OCFD2d_BoundaryConditions.f90 to find the supported BCs"
               stop
            end select
         end

!------------Set parameters --------------------------------
         subroutine set_parameters
            use flow_para
            implicit none

            hx = 1.d0/(nx_global - 1.d0)
            hy = 1.d0/(ny_global - 1.d0)
            hz = 1.d0/(nz_global - 1.d0)

            Cv = 1.d0/(Para%gamma*(Para%gamma - 1.d0)*Para%Ma*Para%Ma)
            Cp = Cv*Para%gamma

!---------Set Scheme%bound_index       (指示6个边界 是否采用边界格式）
! Scheme%Scheme_boundary(:)==-1    ! Ghost-Cell type boundary   (Do not use boundary scheme)
            Scheme%bound_index(:, :) = 0             ! default :  not use boundary scheme

            if (npx .eq. 0 .and. Para%Iperiodic_X .eq. 0) Scheme%bound_index(1, 1) = 1    !  i-
            if (npx .eq. npx0 - 1 .and. Para%Iperiodic_X .eq. 0) Scheme%bound_index(2, 1) = 1         !  i+

            if (npy .eq. 0 .and. Para%Iperiodic_Y .eq. 0) Scheme%bound_index(1, 2) = 1      !  j-
            if (npy .eq. npy0 - 1 .and. Para%Iperiodic_Y .eq. 0) Scheme%bound_index(2, 2) = 1             !  j+

            if (npz .eq. 0 .and. Para%Iperiodic_Z .eq. 0) Scheme%bound_index(1, 3) = 1      !  k-
            if (npz .eq. npz0 - 1 .and. Para%Iperiodic_Z .eq. 0) Scheme%bound_index(2, 3) = 1             !  k+

! ----[Ka1,Kb1]:  Stencil of positive flux F(i+1/2) : [i+Ka1, i+Kb1] ;
!     [Ka2,Kb2]:  Stencil of negative flux ;

            select case (Scheme%Scheme_Invis)
            case (OCFD_Scheme_WENO5)
               Scheme%Ka1 = -2; Scheme%Kb1 = 2      ! Stencil for F(i+1/2) of WENO5+ : [i-2, ..., i+2]
               Scheme%Ka2 = -1; Scheme%Kb2 = 3      ! Stencil for F(i+1/2) of WENO5- : [i-1, ..., i+3]
            case (OCFD_Scheme_WENO7)
               Scheme%Ka1 = -3; Scheme%Kb1 = 3
               Scheme%Ka2 = -2; Scheme%Kb2 = 4
            case (OCFD_Scheme_OMP6)
               Scheme%Ka1 = -3; Scheme%Kb1 = 4
               Scheme%Ka2 = -3; Scheme%Kb2 = 4
            case (OCFD_Scheme_UD7L)                    ! low-dissipative Upwind Difference scheme
               Scheme%Ka1 = -3; Scheme%Kb1 = 4
               Scheme%Ka2 = -3; Scheme%Kb2 = 4
            case (OCFD_Scheme_Hybrid)                    !  Hybrid scheme
               Scheme%Ka1 = -3; Scheme%Kb1 = 4
               Scheme%Ka2 = -3; Scheme%Kb2 = 4
            case (OCFD_Scheme_USER)
!               call Stencil_Scheme_User(Scheme%Ka1, Scheme%Kb1,   Scheme%Ka2, Scheme%Kb2)

            case default
               print *, "The Inviscous Scheme is not supported"
               stop
            end select

!       [Ka,Kb] for hybrid scheme
            Scheme%Ka1_H1 = -3; Scheme%Kb1_H1 = 4; Scheme%Ka2_H1 = -3; Scheme%Kb2_H1 = 4      ! UD7L
            Scheme%Ka1_H2 = -3; Scheme%Kb1_H2 = 3; Scheme%Ka2_H2 = -2; Scheme%Kb2_H2 = 4      ! WENO7
            Scheme%Ka1_H3 = -2; Scheme%Kb1_H3 = 2; Scheme%Ka2_H3 = -1; Scheme%Kb2_H3 = 3      ! WENO5

            call set_scheme_para

         end

!----------------change to capital letters -------------------
         subroutine capitalize(string)
            implicit none
            character(len=*):: string
            integer::i, length
            length = len(string)
            do i = 1, length
               if (lge(string(i:i), 'a') .and. lle(string(i:i), 'z')) then
                  string(i:i) = achar(iachar(string(i:i)) - 32)
               end if
            end do
         end
! ---------Coefficients for low-dissipative upwind scheme (UD7L)
         subroutine set_scheme_para
            use flow_para
            implicit none
            real(kind=OCFD_REAL_KIND):: UD7(8), CD8(8)
            UD7 = (/-3.d0, 25.d0, -101.d0, 319.d0, 214.d0, -38.d0, 4.d0, 0.d0/)       ! coefficients for UD7
            CD8 = (/-3.0, 29.0, -139.0, 533.0, 533.0, -139.0, 29.0, -3.0/)        ! coefficients for CD8
            Scheme%UD7L = Scheme%UD7L_Diss*UD7/420.d0 + (1.d0 - Scheme%UD7L_Diss)*CD8/840.d0
         end

