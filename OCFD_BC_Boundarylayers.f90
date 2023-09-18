! boundary conditions： boundary-layers  (transition/ STBLI/ ...) 
!---------------------------------------------------------------------
   module BC_data 
   use OCFD_precision
   implicit none 
   TYPE bc_type
   integer:: bc_init=0
   Real(kind=OCFD_REAL_KIND),allocatable:: flow_inlet1d(:,:), flow_upper1d(:,:), &
                         flow_inlet2d(:,:,:), flow_upper2d(:,:,:) 
   Real(kind=OCFD_REAL_KIND),allocatable:: wall_pertb(:,:)      ! Wall blow-and-suction
   integer,allocatable:: NonReflect_upper(:,:)          
   
   end TYPE 
   TYPE (bc_type):: BC 
   end    
   
	
	subroutine ocfd_bc_BoundaryLayer
    use flow_data
	use BC_data  
	implicit none 
	integer i,j,k,m
	integer:: BcData_inlet, BcData_upper,bc_upper_nonref,bc_outlet,bc_dis_type,bc_dis_mt,bc_dis_mz
	Real(kind=OCFD_REAL_KIND):: Tw,Wall_Xinit,bc_dis_A,bc_dis_Xbegin,bc_dis_Xend,bc_dis_ZL,bc_dis_freq
	Real(kind=OCFD_REAL_KIND):: u0,v0, ht, tmp,A0
	
!---------------------------------------------------
     BcData_inlet=nint(Para%BC_para(1))      ! inlet: (0 free-stream;  1 1d data;  2 2d data)
	 BcData_upper=nint(Para%Bc_para(2))      ! upper boundary: (0 free-stream, 1 1d data ,  2 2d data)
	 
	 bc_upper_nonref=nint(Para%Bc_para(3))   ! upper boundary: 0 auto,  1 Non-Reflection,  2 Dirichlet 
	 bc_outlet=nint(Para%Bc_para(4))         ! outlet: 0 Non-Reflection,   1 1st order extrapolation,  2 2nd order extrapolation
     Tw=Para%Bc_para(5)                      ! wall temperature
	 Wall_Xinit=Para%Bc_para(6)              ! x location of the wall leading 

! wall perturbation--------	 
	 bc_dis_type=nint(Para%Bc_para(7))     ! wall disturbance type (0 none ;  1 multi-wave blow-and-suction, Ref: Rai MM, AIAA 95-0583)
	 bc_dis_A=Para%Bc_para(8)              ! Amplitude of wall disturbance     
	 bc_dis_Xbegin=Para%Bc_para(9)	       ! Initial location of wall disturbance
	 bc_dis_Xend=Para%Bc_para(10)          ! End location of wall disturbance
	 bc_dis_mt=nint(Para%Bc_para(11))     ! multi-frequency
	 bc_dis_mz=nint(Para%Bc_para(12))     ! multi-wavenumber
	 bc_dis_ZL=Para%Bc_para(13)           ! Spanwise Length
     bc_dis_freq=Para%Bc_para(14)         ! base frequancy of disturbance  
	 
	 
 !-------- j=1 wall ----------------		 
     if(BC%bc_init ==0 ) then 
	  call init_bcdata_boundarylayer (BcData_inlet, BcData_upper)  ! read inlet & upper data 
	  call init_bc_nonRef (bc_upper_nonref)             ! upper boundary:  Non-Reflection or Dirichlet 
      call init_bc_wall_perturbation (bc_dis_type, bc_dis_A, bc_dis_Xbegin, bc_dis_Xend, bc_dis_mt, bc_dis_mz, bc_dis_ZL) 	   
      Bc%bc_init=1                     ! Run only for initial time 
	 endif 
	 
!----------------Inlet boundary -------------------------------------
    u0=cos(Para%AoA) 
    v0=sin(Para%AoA)
  
 if(npx .eq. 0) then 
  if( BcData_inlet == 0) then 
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=1.d0 
    u(1,j,k)=u0
    v(1,j,k)=v0
    w(1,j,k)=0.d0 
    T(1,j,k)=0.d0 
   enddo 
   enddo 
  else if(BcData_inlet == 1) then   
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=BC%flow_inlet1d(j,1) 
    u(1,j,k)=BC%flow_inlet1d(j,2)  
    v(1,j,k)=BC%flow_inlet1d(j,3) 
    w(1,j,k)=BC%flow_inlet1d(j,4) 
    T(1,j,k)=BC%flow_inlet1d(j,5)  
   enddo 
   enddo  	 
  else if (BcData_inlet == 2) then	 
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=BC%flow_inlet2d(j,k,1) 
    u(1,j,k)=BC%flow_inlet2d(j,k,2)  
    v(1,j,k)=BC%flow_inlet2d(j,k,3) 
    w(1,j,k)=BC%flow_inlet2d(j,k,4) 
    T(1,j,k)=BC%flow_inlet2d(j,k,5)  
   enddo 
   enddo  	   
  else 
   print*, "The inlet data is not supported !"
  endif 
 endif 
 
!----------------upper boundary ------------------------------------- 
  if(npy .eq. npy0-1) then 
   do k=1,nz 
   do i=1,nx 
   if(Bc%NonReflect_upper(i,k) .eq. 0) then   ! Dirichlet BC 
   if(BcData_upper .eq. 0) then 
    d(i,ny,k)=1.d0 
    u(i,ny,k)=u0
    v(i,ny,k)=v0
    w(i,ny,k)=0.d0 
    T(i,ny,k)=1.d0 
   else if (BcData_upper .eq. 1) then
    d(i,ny,k)=BC%flow_upper1d(i,1) 
    u(i,ny,k)=BC%flow_upper1d(i,2)    
	v(i,ny,k)=BC%flow_upper1d(i,3) 
    w(i,ny,k)=BC%flow_upper1d(i,4)    
    T(i,ny,k)=BC%flow_upper1d(i,5)
   else if(BcData_upper .eq. 2) then
    d(i,ny,k)=BC%flow_upper2d(i,k,1) 
    u(i,ny,k)=BC%flow_upper2d(i,k,2)    
	v(i,ny,k)=BC%flow_upper2d(i,k,3) 
    w(i,ny,k)=BC%flow_upper2d(i,k,4)    
    T(i,ny,k)=BC%flow_upper2d(i,k,5)  
   endif 
   endif 
   enddo 
   enddo
  endif 
!----wall----------------------------------------------- 
     call get_ht_multifrequancy(ht,tt,bc_dis_mt, bc_dis_freq)
     if(npy.eq.0) then
      do k=1,nz
      do i=1,nx
	   tmp=1.d0/(sqrt(Aix(i,1,k)**2+Aiy(i,1,k)**2+Aiz(i,1,k)**2)) 
	   A0=bc_dis_A*ht* BC%wall_pertb(i,k)
	   u(i,1,k)=A0*Aix(i,1,k)*tmp      !  Aix*tmp  normal vector 
	   v(i,1,k)=A0*Aiy(i,1,k)*tmp
	   w(i,1,k)=A0*Aiz(i,1,k)*tmp
       if(Tw.gt.0) then
        T(i,1,k)=Tw
        d(i,1,k)=(4.d0*d(i,2,k)*T(i,2,k)-d(i,3,k)*T(i,3,k))/(3.d0*Tw)    ! P1=(4.d0*P2-P3)/3.d0  2nd order 
       else
        T(i,1,k)=(4.d0*T(i,2,k)-T(i,3,k))/3.d0
        d(i,1,k)=(4.d0*d(i,2,k)-d(i,3,k))/3.d0
       endif
	  enddo 
	  enddo 
	  
	 endif 
!-------outlet bounary ------------------------------------------
   if(npx .eq. npx0-1) then 
   if(bc_outlet==1) then  
	 do k=1,nz 
	 do j=1,ny 
	  d(nx,j,k)=d(nx-1,j,k)
	  u(nx,j,k)=u(nx-1,j,k)
	  v(nx,j,k)=v(nx-1,j,k)
	  w(nx,j,k)=w(nx-1,j,k)
	  T(nx,j,k)=T(nx-1,j,k)
	 enddo 
	 enddo 
    else if (bc_outlet==2) then 
	 do k=1,nz 
	 do j=1,ny 
	  d(nx,j,k)=2.d0*d(nx-1,j,k)-d(nx-2,j,k)
	  u(nx,j,k)=2.d0*u(nx-1,j,k)-u(nx-2,j,k)
	  v(nx,j,k)=2.d0*v(nx-1,j,k)-v(nx-2,j,k)
	  w(nx,j,k)=2.d0*w(nx-1,j,k)-w(nx-2,j,k)
	  T(nx,j,k)=2.d0*T(nx-1,j,k)-T(nx-2,j,k)
	 enddo 
	 enddo 	
   endif 	
   endif 

    end 	
	
!---------------------------------------------------------------------	
! read  boundary data (inlet & upper) 
 subroutine init_bcdata_boundarylayer (BcData_inlet, BcData_upper)
  use flow_data
  use BC_data 
  implicit none  
  integer::  BcData_inlet, BcData_upper
  integer::  i,j,k,m,i1,j1,k1,ierr
  Real(kind=OCFD_REAL_KIND):: tmp 
  Real(kind=OCFD_REAL_KIND),allocatable:: flow1d0(:,:),flow2d0(:,:)

  
!  BcData_inlet=nint(Para%BC_para(1))      ! inlet type (1 free-stream;  2 1d data;  2 2d data)
!  BcData_upper=nint(Para%Bc_para(2))
  


!----------------read inlet data file  (1d formatted or 2d unformatted )  
 if(BcData_inlet == 1) then        ! 1d inlet data file  
   allocate(BC%flow_inlet1d(ny,5), flow1d0(ny_global,4))
   if(my_id .eq. 0) then 
	open(99,file="flow1d-inlet.dat")
	read(99,*) 
	do j=1,ny_global
	read(99,*) tmp, (flow1d0(j,k),k=1,4)        ! d,u,v,T
	enddo 
    close(99)
    print*, "read 1d inlet data OK"
   endif 
   
   call MPI_bcast(flow1d0,ny_global*4,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)   
   do j=1,ny 
   j1=j_offset(npy)+j-1
   BC%flow_inlet1d(j,1:3)=flow1d0(j1,1:3) 
   BC%flow_inlet1d(j,4)=0.d0     ! w 
   BC%flow_inlet1d(j,5)=flow1d0(j1,4)
   enddo 
   deallocate(flow1d0) 
  
  else if(BcData_inlet == 2)  then   !2d inlet data file 
   allocate(BC%flow_inlet2d(ny,nz,5), flow2d0(ny_global,nz_global))
   if(my_id .eq. 0)  open(99,file="flow2d-inlet.dat",form="unformatted")
    do m=1,5
    if(my_id .eq. 0) read(99) flow2d0 
	
    call MPI_bcast(flow2d0,ny_global*nz_global ,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 	
	do k=1,nz 
	do j=1,ny 
	 j1=j_offset(npy)+j-1
	 k1=k_offset(npz)+k-1 
	 BC%flow_inlet2d(j,k,m)=flow2d0(j1,k1) 
	enddo 
	enddo 
	enddo 
	if(my_id .eq. 0) then 
	  close(99)
	  print*, "read 2d inlet data OK"
	endif 
    deallocate(flow2d0) 
  endif 
  
 !----------------read upper-boundary data file  (1d formatted or 2d unformatted )  
  if(BcData_upper == 1) then        ! 1d upper-boundary data file  
    allocate(BC%flow_upper1d(nx,5), flow1d0(nx_global,4))
   if(my_id .eq. 0) then 
    open(99,file="flow1d-upper.dat")
	read(99,*) 
	do i=1,nx_global
	read(99,*) tmp, (flow1d0(i,k),k=1,4)        ! d,u,v,T
	enddo 
    close(99)
	print*, "read 1d upper boundary data OK"
   endif 
   
   call MPI_bcast(flow1d0,nx_global*4,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)   
   do i=1,nx 
   i1=i_offset(npx)+i-1
   BC%flow_upper1d(i,1:3)=flow1d0(i1,1:3) 
   BC%flow_upper1d(i,4)=0.d0     ! w 
   BC%flow_upper1d(i,5)=flow1d0(i1,4)
   enddo 
   deallocate(flow1d0) 
  
  else if(BcData_upper == 2)  then   !2d boundary data file 
   allocate(BC%flow_upper2d(nx,nz,5), flow2d0(nx_global,nz_global))
    if(my_id .eq. 0)  open(99,file="flow2d-upper.dat",form="unformatted")
    do m=1,5
    if(my_id .eq. 0) read(99) flow2d0 
	
    call MPI_bcast(flow2d0,nx_global*nz_global ,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 	
	do k=1,nz 
	do i=1,nx 
	 i1=i_offset(npx)+i-1
	 k1=k_offset(npz)+k-1 
	 BC%flow_upper2d(i,k,m)=flow2d0(i1,k1) 
	enddo 
	enddo 
	enddo 
	if(my_id .eq. 0) then 
	 close(99)
     print*, "read 2d upper boundary data OK"
	endif
	deallocate(flow2d0) 
  endif  
 end 
 
 
 !-----------------------------------------------------
  subroutine init_bc_nonRef (Bc_NonReflect_upper)
  use flow_data
  use BC_data 
  implicit none  
  integer::  Bc_NonReflect_upper,i,j,k
  Real(kind=OCFD_REAL_KIND)::  sy ,sy0
  Real(kind=OCFD_REAL_KIND),parameter:: epsl=1.d-6 
  
!  Bc_NonReflect_upper=nint(Para%Bc_para(3))      ! 0 Auto,  1 Non-Reflection,  2 Dirichlet boudary
  
  sy0=tan(Para%AoA)
  allocate(BC%NonReflect_upper(nx,nz)) 
  Bc%NonReflect_upper=0
  
  
!  if(npy .eq. 0) then                           ! Bug removed 2021-5-4
  if(npy .eq. npy0-1) then 
    if(Bc_NonReflect_upper == 1) then 
	  BC%NonReflect_upper=1            ! Non-Reflection boundary 
   else if (Bc_NonReflect_upper == 2) then  
      BC%NonReflect_upper=0            ! Dirichlet boundary 
   else if (Bc_NonReflect_upper == 0) then      ! Auto type (Non-Reflection /Dirichlet)
    do k=1,nz 
	do i=1,nx 
	  if(i==1 .and. npx .eq. 0) then
	   sy=(Ayy(i+1,ny,k)-Ayy(i,ny,k))/(Axx(i+1,ny,k)-Axx(i,ny,k))
	  else 
	   sy=(Ayy(i,ny,k)-Ayy(i-1,ny,k))/(Axx(i,ny,k)-Axx(i-1,ny,k))
      endif 
	 if(sy >= sy0+ epsl) then 
	   Bc%NonReflect_upper(i,k)=0       ! Dirichlet BC 
	 else  
	  Bc%NonReflect_upper(i,k)=1        ! Non-Reflection BC
	 endif 
	enddo
	enddo 
   endif    
  endif 
  
end 	
	
   subroutine init_bc_wall_perturbation (bc_dis_type, bc_dis_A, bc_dis_Xbegin, bc_dis_Xend, bc_dis_mt, bc_dis_mz, bc_dis_ZL) 
! TM(k), Amplitude; faiz,fait: phase (random); mzmax, mtmax: wavenumbers
! consider wall jet
   use flow_data
   use BC_data
   implicit none  
   integer i,j,k,m,ierr
   real(kind=OCFD_REAL_KIND),allocatable :: faiz(:),Zl(:)
   real(kind=OCFD_REAL_KIND),parameter:: PI=3.14159265358979d0
   real(kind=OCFD_REAL_KIND):: ztmp,seta,fx,gz,rtmp
   integer:: bc_dis_type,bc_dis_mt,bc_dis_mz
   real(kind=OCFD_REAL_KIND):: bc_dis_A,bc_dis_Xbegin,bc_dis_Xend,bc_dis_ZL
   
   
!	 bc_dis_type=nint(Para%Bc_para(7))     ! wall disturbance type (0 none ;  1 multi-wave blow-and-suction, Ref: Rai MM, AIAA 95-0583)
!	 bc_dis_A=Para%Bc_para(8)              ! Amplitude of wall disturbance     
!	 bc_dis_Xbegin=Para%Bc_para(9)	       ! Initial location of wall disturbance
!	 bc_dis_Xend=Para%Bc_para(10)          ! End location of wall disturbance
!	 bc_dis_mt=nint(Para%Bc_para(11))     ! multi-frequency
!	 bc_dis_mz=nint(Para%Bc_para(12))     ! multi-wavenumber
!	 bc_dis_ZL=Para%Bc_para(13)           ! Spanwise Length
 
	 allocate(BC%wall_pertb(nx,nz))
     BC%wall_pertb=0.d0 
	 
   if(bc_dis_type == 0) then 
      BC%wall_pertb=0.d0       ! no disturbation
 !-----------------------------------------------
   else if (	bc_dis_type == 1) then 	 ! multi-wavenumber perturbation (Rai MM AIAA 95-0583) 
	 
	 if(bc_dis_mz >0) then  
	  allocate(faiz(bc_dis_mz),Zl(bc_dis_mz))
      ztmp=0.d0
      do k=1,bc_dis_mz
      call random_number(faiz(k))
      if(k.eq.1) then
       Zl(k)=1.d0
      else
       zl(k)=zl(k-1)/1.25d0
      endif
       ztmp=ztmp+Zl(k)
      enddo
      do k=1,bc_dis_mz
       zl(k)=zl(k)/ztmp
      enddo
      call MPI_bcast(faiz(1),bc_dis_mz,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     endif 

    do k=1,nz
	do i=1,nx 
	 if(Axx(i,1,k) >=bc_dis_Xbegin .and. Axx(i,1,k) <= bc_dis_Xend) then 
	   seta=2.d0*PI*(Axx(i,1,k)-bc_dis_Xbegin)/(bc_dis_Xend-bc_dis_Xbegin)
	   fx=4.d0/sqrt(27.d0)*sin(seta)*(1.d0-cos(seta))
	 else 
       fx=0.d0 
     endif 
   
     gz=0.d0
     seta=Azz(i,1,k)/bc_dis_ZL
     if(bc_dis_mz > 0) then
	  do m=1,bc_dis_mz
       gz=gz+Zl(m)*sin(2.d0*PI*m* (seta+faiz(m)) )
      enddo
     else if(bc_dis_mz == 0) then
      gz=1.d0
	 else if(bc_dis_mz < 0) then
      gz=sin(-2.d0*PI*bc_dis_mz*seta)
	 endif
	 BC%wall_pertb(i,k)=fx*gz
	enddo
    enddo 
    if(bc_dis_mz >0)  deallocate(faiz,Zl)
!----------------------------------------------- 
  else if (	bc_dis_type == 2) then 	       ! random disturbance 
	do k=1,nz 
	do i=1,nx 
	 if(Axx(i,1,k) >=bc_dis_Xbegin .and. Axx(i,1,k) <= bc_dis_Xend) then         ! a Bug removed (2021-6-21)
	  call random_number(rtmp)
	  BC%wall_pertb(i,k)=2.d0*(rtmp-0.5d0)          ! -1< wall_pertb < 1	  
	 else 
	  BC%wall_pertb(i,k)=0.d0 
	 endif  
    enddo 
	enddo 
  endif 
  end  

  
!--------------------------------	
! perturbation of Rai, see:  Rai MM, AIAA 95-0583
     subroutine get_ht_multifrequancy(ht,tt,mtmax,beta)
     use OCFD_constants
	 use Para_mpi
     implicit none
     integer mtmax,k,m,ierr
     integer,save:: Kflag=0
     real(kind=OCFD_REAL_KIND):: ht,tt,beta,Ttmp,rand_x
     real(kind=OCFD_REAL_KIND),parameter:: PI=3.14159265358979d0
     real(kind=OCFD_REAL_KIND),allocatable,save:: TM(:),fait(:)          ! Amplitute and random phase angle
     
     if(Kflag .eq. 0) then
 	 Kflag=1
     allocate(TM(mtmax),fait(mtmax))

	 Ttmp=0.d0
     do k=1,mtmax
      call random_number(rand_x)
      fait(k)=rand_x
     if(k.eq.1) then
      TM(k)=1.d0
     else
      TM(k)=TM(k-1)/1.25d0
     endif
      Ttmp=Ttmp+TM(k)
     enddo
      do k=1,mtmax
        TM(k)=TM(k)/Ttmp
      enddo
     call MPI_bcast(fait(1),mtmax,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     endif

     ht=0.d0
     if(beta .gt. 0.d0) then      
      do m=1,mtmax
       ht=ht+TM(m)*sin(m*beta*tt+2.d0*PI*fait(m))
       enddo
      else
       ht=1.d0
      endif
      end
	  
!------------modified 2021-5-10 ---------------------
	subroutine ocfd_bc_swept_corner
    use flow_data
	use BC_data  
	implicit none 
	integer i,j,k,m
	integer:: BcData_inlet, BcData_upper,bc_upper_nonref,bc_outlet,bc_dis_type,bc_dis_mt,bc_dis_mz, &
	         bc_Z1, bc_Z2, bc_Wallext
	Real(kind=OCFD_REAL_KIND):: Tw,Wall_Xinit,bc_dis_A,bc_dis_Xbegin,bc_dis_Xend,bc_dis_ZL,bc_dis_freq, &
	                            bc_Wall_ZL
	Real(kind=OCFD_REAL_KIND):: u0,v0, ht, tmp,A0
	Real(kind=OCFD_REAL_KIND),parameter :: epsl=1.d-7
!-----------------------------------------------------------------
     BcData_inlet=1        ! inlet: (0 free-stream;  1 1d data;  2 2d data)
	 BcData_upper=0      ! upper boundary: (0 free-stream, 1 1d data ,  2 2d data)
	 bc_upper_nonref=0   ! upper boundary: 0 auto,  1 Non-Reflection,  2 Dirichlet 

	 bc_outlet=nint(Para%Bc_para(1))         ! outlet: 0 Non-Reflection,   1 1st order extrapolation,  2 2nd order extrapolation
     Tw=Para%Bc_para(2)                      ! wall temperature
     bc_Wall_ZL=Para%Bc_para(3)         ! Spanwise Length of  wall
	 bc_dis_ZL=bc_Wall_ZL
! wall perturbation--------	 
	 bc_dis_type=nint(Para%Bc_para(4))     ! wall disturbance type (0 none ;  1 multi-wave blow-and-suction, Ref: Rai MM, AIAA 95-0583)
	 bc_dis_A=Para%Bc_para(5)              ! Amplitude of wall disturbance     
	 bc_dis_Xbegin=Para%Bc_para(6)	       ! Initial location of wall disturbance
	 bc_dis_Xend=Para%Bc_para(7)          ! End location of wall disturbance
	 bc_dis_mt=nint(Para%Bc_para(8))     ! multi-frequency
	 bc_dis_mz=nint(Para%Bc_para(9))     ! multi-wavenumber
     bc_dis_freq=Para%Bc_para(10)         ! base frequancy of disturbance  
	 bc_Z1=nint(Para%Bc_para(11))  ! Lift side (Z-) boundary; (0 slide wall;  1 exterploation)  
	 bc_Z2=nint(Para%Bc_para(12))  ! Right side (Z+) boundary
!	 bc_Wallext=nint(Para%Bc_para(13))  ! Wall extent region (j=1;  z<0 or z> Zwall)
 !-------- j=1 wall ----------------		 
     if(BC%bc_init ==0 ) then 
	  call init_bcdata_boundarylayer (BcData_inlet, BcData_upper)  ! read inlet & upper data 
	  call init_bc_nonRef (bc_upper_nonref)             ! upper boundary:  Non-Reflection or Dirichlet 
      call init_bc_wall_perturbation (bc_dis_type, bc_dis_A, bc_dis_Xbegin, bc_dis_Xend, bc_dis_mt, bc_dis_mz, bc_dis_ZL) 	   
      Bc%bc_init=1                     ! Run only for initial time 
	 endif 
	 
!----------------Inlet boundary -------------------------------------
    u0=cos(Para%AoA) 
    v0=sin(Para%AoA)
  
  if(npx .eq. 0) then 
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=BC%flow_inlet1d(j,1) 
    u(1,j,k)=BC%flow_inlet1d(j,2)  
    v(1,j,k)=BC%flow_inlet1d(j,3) 
    w(1,j,k)=BC%flow_inlet1d(j,4) 
    T(1,j,k)=BC%flow_inlet1d(j,5)  
   enddo 
   enddo  
  endif 
 
!----------------upper boundary ------------------------------------- 
  if(npy .eq. npy0-1) then 
   do k=1,nz 
   do i=1,nx 
   if(Bc%NonReflect_upper(i,k) .eq. 0) then   ! Dirichlet BC 
    d(i,ny,k)=1.d0 
    u(i,ny,k)=u0
    v(i,ny,k)=v0
    w(i,ny,k)=0.d0 
    T(i,ny,k)=1.d0 
    endif 
   enddo 
   enddo
  endif 
!----wall----------------------------------------------- 
     call get_ht_multifrequancy(ht,tt,bc_dis_mt, bc_dis_freq)
     if(npy.eq.0) then
      do k=1,nz
      do i=1,nx
      if( (Azz(i,1,k) >= -epsl .and. Azz(i,1,k) <=  bc_Wall_ZL+epsl) .or. Ayy(i,1,k) <= epsl ) then   ! non-slide wall , modifyied 2021-5-10	  
	   tmp=1.d0/(sqrt(Aix(i,1,k)**2+Aiy(i,1,k)**2+Aiz(i,1,k)**2)) 
	   A0=bc_dis_A*ht* BC%wall_pertb(i,k)
	   u(i,1,k)=A0*Aix(i,1,k)*tmp      !  Aix*tmp  normal vector 
	   v(i,1,k)=A0*Aiy(i,1,k)*tmp
	   w(i,1,k)=A0*Aiz(i,1,k)*tmp
       if(Tw.gt.0) then
        T(i,1,k)=Tw
        d(i,1,k)=(4.d0*d(i,2,k)*T(i,2,k)-d(i,3,k)*T(i,3,k))/(3.d0*Tw)    ! P1=(4.d0*P2-P3)/3.d0  2nd order 
       else
        T(i,1,k)=(4.d0*T(i,2,k)-T(i,3,k))/3.d0
        d(i,1,k)=(4.d0*d(i,2,k)-d(i,3,k))/3.d0
       endif
	  else    
    ! 1st order exterpolation   
 	   d(i,1,k)=d(i,2,k)
	   T(i,1,k)=T(i,2,k)
	   u(i,1,k)=u(i,2,k)
       v(i,1,k)=v(i,2,k)        ! modified 2021-4-26 
	   w(i,1,k)=w(i,2,k) 
	  endif
	  enddo 
	  enddo 
	 endif 
	 
!-------outlet bounary ------------------------------------------
   if(npx .eq. npx0-1) then 
   if(bc_outlet==1) then  
	 do k=1,nz 
	 do j=1,ny 
	  d(nx,j,k)=d(nx-1,j,k)
	  u(nx,j,k)=u(nx-1,j,k)
	  v(nx,j,k)=v(nx-1,j,k)
	  w(nx,j,k)=w(nx-1,j,k)
	  T(nx,j,k)=T(nx-1,j,k)
	 enddo 
	 enddo 
    else if (bc_outlet==2) then 
	 do k=1,nz 
	 do j=1,ny 
	  d(nx,j,k)=2.d0*d(nx-1,j,k)-d(nx-2,j,k)
	  u(nx,j,k)=2.d0*u(nx-1,j,k)-u(nx-2,j,k)
	  v(nx,j,k)=2.d0*v(nx-1,j,k)-v(nx-2,j,k)
	  w(nx,j,k)=2.d0*w(nx-1,j,k)-w(nx-2,j,k)
	  T(nx,j,k)=2.d0*T(nx-1,j,k)-T(nx-2,j,k)
	 enddo 
	 enddo 	
   endif 	
   endif 

!----------side boundary  (slide (symmetry)/exterpolation )-------------------------
   if(npz .eq. 0 ) then 
!   if(bc_Z1 .eq. 0) then   ! not support 2021-6-14  (do not work well)
!	do j=1,ny 
!    do i=1,nx 
!    d(i,j,1)=d(i,j,2) 
!    u(i,j,1)=u(i,j,2)
!    v(i,j,1)=v(i,j,2)
!    w(i,j,1)=0.d0 
!    T(i,j,1)=T(i,j,2)
!    enddo 
!    enddo  
!    else if (bc_Z1 .eq. 1) then  
  if (bc_Z1 .eq. 1) then 
	do j=1,ny 
    do i=1,nx 
     d(i,j,1)=d(i,j,2) 
     u(i,j,1)=u(i,j,2)
     v(i,j,1)=v(i,j,2)
     w(i,j,1)=w(i,j,2)
     T(i,j,1)=T(i,j,2)
    enddo 
    enddo 
   endif 
   endif 
   
   if(npz .eq. npz0-1 ) then 
!   if(bc_Z2 .eq. 0) then   
!   do j=1,ny 
!   do i=1,nx 
!    d(i,j,nz)=d(i,j,nz-1) 
!    u(i,j,nz)=u(i,j,nz-1)
!    v(i,j,nz)=v(i,j,nz-1)
!    w(i,j,nz)=0.d0
!    T(i,j,nz)=T(i,j,nz-1)
!   enddo   
!   enddo 
!   else if(bc_Z2 .eq. 1) then
 if(bc_Z2 .eq. 1) then
   do j=1,ny 
   do i=1,nx 
    d(i,j,nz)=d(i,j,nz-1) 
    u(i,j,nz)=u(i,j,nz-1)
    v(i,j,nz)=v(i,j,nz-1)
    w(i,j,nz)=w(i,j,nz-1)
    T(i,j,nz)=T(i,j,nz-1)
   enddo   
   enddo
   endif 	
   endif   
   end 	
	
	  
	  