! -------------------------------------------------------------------------------------------------
! OpenCFD-SC Ver 2.x,   3D Compressible Naver-Stoker Finite Difference solver
! Copyright by Li Xinliang, LHD, Institute of Mechanics.   Email: lixl@imech.ac.cn
! Last version 2021-2
! New feature in version 2:
!   Character flux resconstruction is supported
!   WENO5-type boundary scheme
!   code is more simple (3d Jocabian solver is used for all geometry
!------------------------------------------------------------------------------------------------------
!  OpenCFD-SC use Double-Precision for default, If you want SINGLE-PRECISION computation,
!  you can change OCFD_REAL_KIND=8 to OCFD_REAL_KIND=4, and change
!  OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION to OCFD_DATA_TYPE=MPI_REAL in OCFD_parameters.f90
!-------------------------------------------------------------------------------------------------------
  program opencfd_sc
     use flow_para
     implicit none
     integer, parameter ::IBUFFER_SIZE = 10000000
     real(kind=8):: BUFFER_MPI(IBUFFER_SIZE)
     integer:: ierr
!----------------------------------------------
     call mpi_init(ierr)                       ! initial of MPI
     call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierr)   ! get my_id
     call MPI_BUFFER_ATTACH(BUFFER_MPI, 8*IBUFFER_SIZE, ierr)    ! attach buffer for MPI_Bsend( )
!----------------------------

     call read_parameter      ! read opencfd2.in (flow parameters)
     call partation3d_mpi     ! partation for MPI  (npx0*npy0*npz0  MPI precesses)
     call set_parameters      ! Set parameters such as Scheme%Bound_Index
     call NS_Solver           ! Solve 3d N-S eq.
     call mpi_finalize(ierr)
  end program opencfd_sc

