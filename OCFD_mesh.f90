  subroutine init_mesh
  Use flow_data 
  implicit none
  integer:: i,j,k
  
    if(Para%Iflag_Gridtype .eq. GRID1D) then 
	 call read_mesh1d
	else if(Para%Iflag_Gridtype .eq. GRID2D_PLANE) then
     call read_mesh2d_plane
    else if(Para%Iflag_Gridtype .eq. GRID2D_AXIAL_SYMM) then 
	 call read_mesh2d_AxialSymm
  	else if(Para%Iflag_Gridtype .eq. GRID3D .or. Para%Iflag_Gridtype .eq. GRID_AND_JACOBIAN3D ) then
	 call read_mesh3d 
   	endif 
         
	   call exchange_boundary_xyz(Axx)
       call exchange_boundary_xyz(Ayy)
       call exchange_boundary_xyz(Azz)
    
       call updata_Axyz_periodic
	   
	   if(Para%Iflag_Gridtype  .ne. GRID_AND_JACOBIAN3D ) then
	     call comput_Jacobian3d
       endif 
	   
 
       call exchange_boundary_xyz(Akx)
       call exchange_boundary_xyz(Aky)
       call exchange_boundary_xyz(Akz)
       call exchange_boundary_xyz(Aix)
       call exchange_boundary_xyz(Aiy)
       call exchange_boundary_xyz(Aiz)
       call exchange_boundary_xyz(Asx)
       call exchange_boundary_xyz(Asy)
       call exchange_boundary_xyz(Asz)
       call exchange_boundary_xyz(Ajac)
 
       call Jac_Ghost_boundary   ! Ghost Cell boundary for Jacobian coefficients
 
	   do k=1-LAP,nz+LAP
	   do j=1-LAP,ny+LAP
	   do i=1-LAP,nx+LAP
	    Akx1(i,j,k)=Akx(i,j,k)/Ajac(i,j,k) 
	    Aky1(i,j,k)=Aky(i,j,k)/Ajac(i,j,k) 
	    Akz1(i,j,k)=Akz(i,j,k)/Ajac(i,j,k)
	    Aix1(i,j,k)=Aix(i,j,k)/Ajac(i,j,k) 
	    Aiy1(i,j,k)=Aiy(i,j,k)/Ajac(i,j,k) 
	    Aiz1(i,j,k)=Aiz(i,j,k)/Ajac(i,j,k)
	    Asx1(i,j,k)=Asx(i,j,k)/Ajac(i,j,k) 
	    Asy1(i,j,k)=Asy(i,j,k)/Ajac(i,j,k) 
	    Asz1(i,j,k)=Asz(i,j,k)/Ajac(i,j,k)
      enddo 
	  enddo
	  enddo 
      if(my_id ==0) print*, "Initial Mesh OK "

      end
!--------------------------------------------------------
	     	
  subroutine comput_Jacobian3d
   Use flow_data
   implicit none
   real(kind=OCFD_REAL_KIND), allocatable,dimension(:,:,:):: xi,xj,xk, yi, yj, yk, zi, zj, zk
   real(kind=OCFD_REAL_KIND):: xi1,xj1,xk1, yi1, yj1, yk1 , zi1, zj1, zk1, Jac1
   integer:: i,j,k 

	 allocate (xi(nx,ny,nz),xj(nx,ny,nz),xk(nx,ny,nz),   &
               yi(nx,ny,nz),yj(nx,ny,nz),yk(nx,ny,nz),   &
	           zi(nx,ny,nz),zj(nx,ny,nz),zk(nx,ny,nz))

	 call OCFD_dx0(Axx,xi,Scheme%Scheme_Vis)
	 call OCFD_dx0(Ayy,yi,Scheme%Scheme_Vis)
     call OCFD_dx0(Azz,zi,Scheme%Scheme_Vis)
     call OCFD_dy0(Axx,xj,Scheme%Scheme_Vis)
     call OCFD_dy0(Ayy,yj,Scheme%Scheme_Vis)
     call OCFD_dy0(Azz,zj,Scheme%Scheme_Vis)
     call OCFD_dz0(Axx,xk,Scheme%Scheme_Vis)
     call OCFD_dz0(Ayy,yk,Scheme%Scheme_Vis)
     call OCFD_dz0(Azz,zk,Scheme%Scheme_Vis)
  

	 do k=1,nz
	 do j=1,ny
	 do i=1,nx
	  xi1=xi(i,j,k); xj1=xj(i,j,k); xk1=xk(i,j,k)
	  yi1=yi(i,j,k); yj1=yj(i,j,k); yk1=yk(i,j,k)
	  zi1=zi(i,j,k); zj1=zj(i,j,k); zk1=zk(i,j,k)
	  Jac1=1.d0/(xi1*yj1*zk1+yi1*zj1*xk1+zi1*xj1*yk1-zi1*yj1*xk1-yi1*xj1*zk1-xi1*zj1*yk1)   ! 1./Jocabian = d(x,y,z)/d(i,j,k) 
      Ajac(i,j,k)=Jac1
	  Akx(i,j,k)=Jac1*(yj1*zk1-zj1*yk1)
	  Aky(i,j,k)=Jac1*(zj1*xk1-xj1*zk1)
	  Akz(i,j,k)=Jac1*(xj1*yk1-yj1*xk1)
	  Aix(i,j,k)=Jac1*(yk1*zi1-zk1*yi1)
	  Aiy(i,j,k)=Jac1*(zk1*xi1-xk1*zi1)
	  Aiz(i,j,k)=Jac1*(xk1*yi1-yk1*xi1)
	  Asx(i,j,k)=Jac1*(yi1*zj1-zi1*yj1)
	  Asy(i,j,k)=Jac1*(zi1*xj1-xi1*zj1)
	  Asz(i,j,k)=Jac1*(xi1*yj1-yi1*xj1)
      if(Jac1 .lt. 0) then 
	   print*, " Jocabian < 0 !!! , Jac=", Jac1  
	   print*, "i,j,k=", i_offset(npx)+i-1, j_offset(npy)+j-1, k_offset(npz)+k-1
  	  endif 
	enddo
	enddo
	enddo
	deallocate ( xi,xj,xk, yi, yj, yk, zi, zj, zk)
   end

	subroutine updata_Axyz_periodic
	 use flow_data 
	 implicit none 
	 integer:: i,j,k  	 
	 if(Para%Iperiodic_X .eq. 1 .and. npx .eq. 0 ) then 
	  do k=1,nz
	  do j=1,ny 
	  do i=1-LAP,0
	   Axx(i,j,k)=Axx(i,j,k)-Para%Periodic_ISpan(1)
	   Ayy(i,j,k)=Ayy(i,j,k)-Para%Periodic_ISpan(2)
	   Azz(i,j,k)=Azz(i,j,k)-Para%Periodic_ISpan(3)	   
	  enddo
	  enddo 
	  enddo 
	 endif
	 
	 if(Para%Iperiodic_X .eq. 1 .and. npx .eq. npx0-1 ) then 
	  do k=1,nz 
	  do j=1,ny 
	  do i=nx+1,nx+LAP
	   Axx(i,j,k)=Axx(i,j,k)+Para%Periodic_ISpan(1)
	   Ayy(i,j,k)=Ayy(i,j,k)+Para%Periodic_ISpan(2)
	   Azz(i,j,k)=Azz(i,j,k)+Para%Periodic_ISpan(3)	   
	  enddo
	  enddo 
	  enddo 
	 endif	 

	 if(Para%Iperiodic_Y .eq. 1 .and. npy .eq. 0 ) then 
     do k=1,nz 	 
	 do j=1-LAP,0 
	 do i=1,nx
	  Axx(i,j,k)=Axx(i,j,k)-Para%Periodic_JSpan(1)
	  Ayy(i,j,k)=Ayy(i,j,k)-Para%Periodic_JSpan(2)
	  Azz(i,j,k)=Azz(i,j,k)-Para%Periodic_JSpan(3) 
	 enddo
	 enddo 
	 enddo 
	 endif	 
	 
	 if(Para%Iperiodic_Y .eq. 1 .and. npy .eq. npy0-1 ) then 
	  do k=1,nz 
	  do j=ny+1,ny+LAP 
	  do i=1,nx
	   Axx(i,j,k)=Axx(i,j,k)+Para%Periodic_JSpan(1)
	   Ayy(i,j,k)=Ayy(i,j,k)+Para%Periodic_JSpan(2)
	   Azz(i,j,k)=Azz(i,j,k)+Para%Periodic_JSpan(3)
	  enddo
	  enddo 
	  enddo
	 endif	

	 if(Para%Iperiodic_Z .eq. 1 .and. npz .eq. 0 ) then 
     do k=1-LAP,0 	 
	 do j=1,ny 
	 do i=1,nx
	  Axx(i,j,k)=Axx(i,j,k)-Para%Periodic_KSpan(1)
	  Ayy(i,j,k)=Ayy(i,j,k)-Para%Periodic_KSpan(2)
	  Azz(i,j,k)=Azz(i,j,k)-Para%Periodic_KSpan(3) 
	 enddo
	 enddo 
	 enddo 
	 endif	
	 
	 if(Para%Iperiodic_Z .eq. 1 .and. npz .eq. npz0-1 ) then 
	  do k=nz+1,nz+LAP 
	  do j=1,ny 
	  do i=1,nx
	   Axx(i,j,k)=Axx(i,j,k)+Para%Periodic_KSpan(1)
	   Ayy(i,j,k)=Ayy(i,j,k)+Para%Periodic_KSpan(2)
	   Azz(i,j,k)=Azz(i,j,k)+Para%Periodic_KSpan(3)
	  enddo
	  enddo 
	  enddo
	 endif		 
	 
    end 
!--------------------------------------------
	subroutine read_mesh3d
     Use flow_data 
     implicit none	 
      if(my_id .eq. 0) then 
	  print*, "read 3D mesh ..."
      open(56,file='OCFD-grid.dat',form='unformatted')
      endif
	  
      call read_3d1(56,Axx)
      call read_3d1(56,Ayy)
      call read_3d1(56,Azz)
	 
	 if(Para%Iflag_Gridtype .eq. GRID_AND_JACOBIAN3D) then 	 		
      call read_3d1(56,Akx)
      call read_3d1(56,Aky)
      call read_3d1(56,Akz)
      call read_3d1(56,Aix)
      call read_3d1(56,Aiy)
      call read_3d1(56,Aiz)
      call read_3d1(56,Asx)
      call read_3d1(56,Asy)
      call read_3d1(56,Asz)
      call read_3d1(56,Ajac)
     endif 
    
	 if(my_id .eq. 0) then
         print*, 'read 3D mesh OK '
         close(56)
      endif
    end	
!--------------------------------------------
	subroutine read_mesh1d
     Use flow_data 
     implicit none	
     real(kind=OCFD_REAL_KIND), allocatable,dimension(:):: x0, y0, z0 
	 integer:: i, j, k, i1,j1,k1,ierr
     allocate(x0(nx_global),y0(ny_global),z0(nz_global))
	 
      if(my_id .eq. 0) then 
	  print*, "read 1D mesh ..."
      open(56,file='OCFD-grid.dat',form='unformatted')
      read(56) x0 
	  read(56) y0
	  read(56) z0
      close(56)     
	  endif
	  call MPI_bcast(x0(1),nx_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	
	  call MPI_bcast(y0(1),ny_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	  
	  call MPI_bcast(z0(1),nz_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	  
	  do k=1,nz 
	  do j=1,ny 
	  do i=1,nx 
	  i1=i_offset(npx)+i-1
	  j1=j_offset(npy)+j-1
	  k1=k_offset(npz)+k-1
	  Axx(i,j,k)=x0(i1)
      Ayy(i,j,k)=y0(j1)
	  Azz(i,j,k)=z0(k1)
	  enddo 
      enddo
      enddo 
     deallocate(x0,y0,z0)
    end		
!--------------------------------------------
! 2D PLANE type mesh,  2d in x-y plane (nx,ny);  1d in z- direcition (nz)
	subroutine read_mesh2d_plane
     Use flow_data 
     implicit none	
     real(kind=OCFD_REAL_KIND), allocatable:: x2(:,:),y2(:,:),z1(:)
	 integer:: i, j, k, i1,j1,k1,ierr
     allocate(x2(nx_global,ny_global),y2(nx_global,ny_global),z1(nz_global))
	 
      if(my_id .eq. 0) then 
	  print*, "read 1D mesh ..."
      open(56,file='OCFD-grid.dat',form='unformatted')
      read(56) x2 
	  read(56) y2
	  read(56) z1
      close(56)     
	  endif
	  call MPI_bcast(x2(1,1),nx_global*ny_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	
	  call MPI_bcast(y2(1,1),nx_global*ny_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	  
	  call MPI_bcast(z1(1),nz_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	  
	  do k=1,nz 
	  do j=1,ny 
	  do i=1,nx 
	  i1=i_offset(npx)+i-1
	  j1=j_offset(npy)+j-1
	  k1=k_offset(npz)+k-1
	  Axx(i,j,k)=x2(i1,j1)
      Ayy(i,j,k)=y2(i1,j1)
	  Azz(i,j,k)=z1(k1)
	  enddo 
      enddo
      enddo 
     deallocate(x2,y2,z1)
    end		
	
! 2D Axial-Symmetry type mesh,  2d in x-y plane (nx,ny);  1d in z- direcition (nz)
	subroutine read_mesh2d_AxialSymm
     Use flow_data 
     implicit none	
     real(kind=OCFD_REAL_KIND), allocatable:: x2(:,:),R2(:,:),seta1(:)
	 integer:: i, j, k, i1,j1,k1,ierr
     allocate(x2(nx_global,ny_global),R2(nx_global,ny_global),seta1(nz_global))
	 
      if(my_id .eq. 0) then 
	  print*, "read 1D mesh ..."
      open(56,file='OCFD-grid.dat',form='unformatted')
      read(56) x2 
	  read(56) R2
	  read(56) seta1
      close(56)     
	  endif
	  call MPI_bcast(x2(1,1),nx_global*ny_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	
	  call MPI_bcast(R2(1,1),nx_global*ny_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	  
	  call MPI_bcast(seta1(1),nz_global,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	  
	  do k=1,nz 
	  do j=1,ny 
	  do i=1,nx 
	  i1=i_offset(npx)+i-1
	  j1=j_offset(npy)+j-1
	  k1=k_offset(npz)+k-1
	  Axx(i,j,k)=x2(i1,j1)
      Ayy(i,j,k)=R2(i1,j1)*cos(seta1(k))     ! y- is horizontal     ; seta-  angle with y- axil
	  Azz(i,j,k)=R2(i1,j1)*sin(seta1(k))     ! z- is up
	  enddo 
      enddo
      enddo 
     deallocate(x2,R2,seta1)
    end			