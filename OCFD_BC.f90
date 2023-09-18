! boundary condition
     subroutine OCFD_bc
     use flow_data
	 implicit none
	 integer i,j,k
! -------Set boundary condition ----------------------
     
     select case (Para%IBC )
     case (BC_None) 
!       No boundary Condition    (such as periodic boundary in both x- and y- )
     case (BC_BoundaryLayer) 
        call  ocfd_bc_BoundaryLayer
     case (BC_Blunt2d)             ! Double Mach Reflaction flow 
      call ocfd_bc_blunt2d
     
     case (BC_SweptCorner)	 ! Swept Compression Corner
	  call ocfd_bc_swept_corner
	 case (BC_User_Def) 
	  call ocfd_bc_User         ! 用户自定义的边界条件
	 case default 
      if(my_id .eq. 0)  print*, "This boundary condition is not supported!!!"
      stop
     end select 	  

	 
!-----------------------------------------------------------------------------------------	   
!c----updata conversion values in boundary   更新守恒变量的边界值------------------------
     if(npx.eq.0) then
      do k=1,nz
      do j=1,ny
      f(1,j,k,1)=d(1,j,k)
      f(1,j,k,2)=d(1,j,k)*u(1,j,k)
      f(1,j,k,3)=d(1,j,k)*v(1,j,k)
      f(1,j,k,4)=d(1,j,k)*w(1,j,k)
      f(1,j,k,5)=d(1,j,k)*(Cv*T(1,j,k)+(u(1,j,k)**2+v(1,j,k)**2+w(1,j,k)**2)*0.5d0)
      enddo
      enddo
     endif
     if(npx.eq.npx0-1) then
      do k=1,nz
      do j=1,ny
      f(nx,j,k,1)=d(nx,j,k)
      f(nx,j,k,2)=d(nx,j,k)*u(nx,j,k)
      f(nx,j,k,3)=d(nx,j,k)*v(nx,j,k)
      f(nx,j,k,4)=d(nx,j,k)*w(nx,j,k)
      f(nx,j,k,5)=d(nx,j,k)*(Cv*T(nx,j,k)+(u(nx,j,k)**2+v(nx,j,k)**2+w(nx,j,k)**2)*0.5d0)
      enddo
      enddo
    endif
!------------------------------
     if(npy.eq.0) then
      do k=1,nz
      do i=1,nx
      f(i,1,k,1)=d(i,1,k)
      f(i,1,k,2)=d(i,1,k)*u(i,1,k)
      f(i,1,k,3)=d(i,1,k)*v(i,1,k)
      f(i,1,k,4)=d(i,1,k)*w(i,1,k)
      f(i,1,k,5)=d(i,1,k)*(Cv*T(i,1,k)+(u(i,1,k)**2+v(i,1,k)**2+w(i,1,k)**2)*0.5d0)
      enddo
      enddo
     endif

     if(npy.eq.npy0-1) then
      do k=1,nz
      do i=1,nx
      f(i,ny,k,1)=d(i,ny,k)
      f(i,ny,k,2)=d(i,ny,k)*u(i,ny,k)
      f(i,ny,k,3)=d(i,ny,k)*v(i,ny,k)
      f(i,ny,k,4)=d(i,ny,k)*w(i,ny,k)
      f(i,ny,k,5)=d(i,ny,k)*(Cv*T(i,ny,k)+(u(i,ny,k)**2+v(i,ny,k)**2+w(i,ny,k)**2)*0.5d0)
      enddo
      enddo
     endif
!--------------------------
     if(npz.eq.0) then
      do j=1,ny
      do i=1,nx
      f(i,j,1,1)=d(i,j,1)
      f(i,j,1,2)=d(i,j,1)*u(i,j,1)
      f(i,j,1,3)=d(i,j,1)*v(i,j,1)
      f(i,j,1,4)=d(i,j,1)*w(i,j,1)
      f(i,j,1,5)=d(i,j,1)*(Cv*T(i,j,1)+(u(i,j,1)**2+v(i,j,1)**2+w(i,j,1)**2)*0.5d0)
      enddo
      enddo
     endif
     if(npz.eq.npz0-1) then
      do j=1,ny
      do i=1,nx
      f(i,j,nz,1)=d(i,j,nz)
      f(i,j,nz,2)=d(i,j,nz)*u(i,j,nz)
      f(i,j,nz,3)=d(i,j,nz)*v(i,j,nz)
      f(i,j,nz,4)=d(i,j,nz)*w(i,j,nz)
      f(i,j,nz,5)=d(i,j,nz)*(Cv*T(i,j,nz)+(u(i,j,nz)**2+v(i,j,nz)**2+w(i,j,nz)**2)*0.5d0)
      enddo
      enddo
     endif

      call Flow_Ghost_boundary       ! Ghost Cell for flow field (d, u, v, w, T)

  end 	 

