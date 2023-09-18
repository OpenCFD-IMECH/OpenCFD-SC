!-----OpenCFD-SC version 2 --------------------------------------------------------------------
! Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Codey by Li Xinliang, 2021-2
!----------------------------------------------------------------------------------------------

! Double precsion (real*8)  or  Single Precision (real*4)       (default: double precision)
module OCFD_precision
implicit none
include "mpif.h"
!------For Doubleprecision  (real*8)--------------------------------------------------------------------------
integer,parameter::OCFD_REAL_KIND=8,  OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
! ------For Single precision (real*4)-------------------------------------------------------------------------
!     integer,parameter::OCFD_REAL_KIND=4,  OCFD_DATA_TYPE=MPI_REAL             !  single precision computing
 end module OCFD_precision


 !------parameters used in OpenCFD---------------------------------------
module OCFD_constants
 Use OCFD_precision
 implicit none
 integer,parameter:: Nvars=5           ! 5 conservative variables, 5 Equations 
 integer,parameter:: OCFD_Turb_None=0                                              ! Turbulence model
                                
 integer,parameter::  OCFD_Scheme_WENO5=1,OCFD_Scheme_WENO7=2,  &                  ! Schemes 
                      OCFD_Scheme_OMP6=3, OCFD_Scheme_UD7L=4,   &
					  OCFD_Scheme_CD6=5,OCFD_Scheme_CD8=6,   &
					  OCFD_Scheme_Hybrid=10,  &
				      OCFD_Scheme_USER=99
 
 integer,parameter:: OCFD_Split_SW=1, OCFD_Split_LLF=2 
 
 integer,parameter:: BC_None=0,BC_Blunt2d=1,  BC_BoundaryLayer=2, BC_SweptCorner=3, BC_User_Def=99
  
 integer,parameter:: GRID1D=10, GRID2D_PLANE=20, GRID2D_AXIAL_SYMM=21, GRID3D=30, GRID_AND_JACOBIAN3D=31  
 
 integer,parameter:: DEL_LIFT=1, DEL_RIGHT=2, DEL_NONE=0             !   WENO5-type boundary Scheme  
 integer,parameter:: OCFD_ANA_USER=99, OCFD_ANA_time_average=100, &
                     OCFD_ANA_Q=101, OCFD_ANA_BOX=102, OCFD_ANA_SAVEDATA=103, OCFD_ANA_Corner=104
end 
 
!----mpi parameters------------------------------------------------ 
module Para_mpi
 implicit none
	integer:: my_id, npx,npy,npz, npx0,npy0,npz0, np_size, &         ! npx0, npy0 : zone numbers in x- and y- direction
        MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z,   MPI_COMM_XY,MPI_COMM_YZ,  MPI_COMM_XZ, &
		ID_XM1,ID_XP1,ID_YM1,ID_YP1, ID_ZM1, ID_ZP1                  ! ID_XM1:  MPI ID for zone-1 (lift one) in x-direction
    integer:: LAP                                                    ! Overlap length for Block-Block		
 	integer,allocatable,dimension(:):: i_offset,j_offset, k_offset, i_nn,j_nn, k_nn
    integer:: TYPE_LAPX1,TYPE_LAPY1,TYPE_LAPZ1,TYPE_LAPX2,TYPE_LAPY2,TYPE_LAPZ2           ! user defined MPI DATA TYPE (for Send & Recv)
 end 
 
 
!-----------parameters for numerical schemes------------------------
module Scheme_Para
 use OCFD_precision
 implicit none 
  TYPE TYPE_Scheme
   integer:: Scheme_Invis, Scheme_Vis
   integer:: Bound_index(2,3), Scheme_boundary(6)    !Bound_index(:,:)  if using boundary scheme(0/1)
           !Scheme_boundary(:) : boundary scheme type ( i-, i+, j-, j+, k-, k+:  0 : default, WENO5 Del-substencil type; -1 full-WENO5 type)  
   integer:: Ka1,Kb1,Ka2,Kb2                      ! [Ka1,Kb1]  stencil for flux+ ;   [Ka2,Kb2] stencil for flux- 
   integer:: Ka1_H1,Kb1_H1,Ka2_H1,Kb2_H1,Ka1_H2,Kb1_H2,Ka2_H2,Kb2_H2,Ka1_H3,Kb1_H3,Ka2_H3,Kb2_H3       ! [Ka,Kb] for Hybrid scheme (scheme1, scheme2, scheme3)   
   Real(kind=OCFD_REAL_KIND):: UD7L_Diss, UD7L(8) , Hybrid_para(100)        
   ! UD7L_Diss:   dissipation coefficient (0-1;  0 : CD8,   1: UD7)
   ! UD7L coefficients for 7th low-dissipation upwind difference sheme  
  end TYPE 
  TYPE (TYPE_Scheme):: Scheme 

end  


!-----------flow parameters ---------------------------------------------- 
module flow_para
 use OCFD_constants
 use Para_mpi 
 use Scheme_para
 implicit none 
 
 
  TYPE TYPE_Para         ! input parameters 
  integer::   Iperiodic_X,Iperiodic_Y,Iperiodic_Z          ! 1/0 : If periodic or not
  integer::   IF_Scheme_Character , IF_Viscous ,Istep_Show, Istep_Save, Flux_Splitting,&
              Iflag_Gridtype ,  &                          ! 0: Grid only, 1: Grid and Jocabian  
			  IBC , IF_Mass_Force, & ! boundary conditon ;  IF_Mase_Force :  0 (default): no mass force;   1 with mass force  
              ANA_Number, &          ! number of analysis process			  
			  Nfiltering, &          ! Number of filtering
              Ghost_cell(6)          ! 0 : non-GhostCell (default);  1 Ghost Cell (Auto, exterpolation);  2 Ghost cell (user defined)
  Real(kind=OCFD_REAL_KIND):: Periodic_ISpan(3), Periodic_JSpan(3) , Periodic_KSpan(3)		    ! Span in peridoic direction
  Real(kind=OCFD_REAL_KIND):: Re, Ma, gamma, Pr , Ref_Amu_T0 ,dt, End_time, AoA           ! AoA angle of attack
  Real(kind=OCFD_REAL_KIND):: BC_para(100)                    ! parameters for boundary conditions 
  Real(kind=OCFD_REAL_KIND):: Mass_Force(3)                   ! Mass force  (such as gravity)
  Real(kind=OCFD_REAL_KIND):: ANA_Para(100,10)                ! parameters for post analysis code 
  Real(kind=OCFD_REAL_KIND):: Filter(100,10)                ! parameters for filtering 
  end TYPE 
 
 TYPE (TYPE_Para):: Para
 
 integer:: nx_global, ny_global,nz_global,  nx, ny, nz, Istep            
 Real(kind=OCFD_REAL_KIND):: Cp, Cv, tt, hx, hy, hz
end 

!---------flow data---------------------------------------------------------
module flow_data
 use flow_para 
 implicit none
 
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: Axx,Ayy,Azz,  &   ! Coordinates 
           Akx,Aky,Akz,Aix,Aiy,Aiz,Asx,Asy,Asz, Ajac , &       !  Jacobian coefficients
           Akx1,Aky1,Akz1,Aix1,Aiy1,Aiz1,Asx1,Asy1,Asz1          ! Akx1=Akx/Ajac 
 real(kind=OCFD_REAL_KIND),allocatable:: f(:,:,:,:),fn(:,:,:,:), du(:,:,:,:), Amu(:,:,:),Amu_t(:,:,:)          
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: d,u,v,w,T
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: Rhybrid         ! Index of hybrid scheme  
end 