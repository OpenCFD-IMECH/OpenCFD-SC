! boundary conditions defined 已预设的几种边界条件
	
!--------Blunt wedge (or half-cylinder )--------------	
	subroutine ocfd_bc_blunt2d
    use flow_data
	implicit none 
	integer i,j,k
!-------- j=1 wall ----------------	
     if(npy.eq.0) then
      do k=1,nz
      do i=1,nx
	   d(i,1,k)=d(i,2,k)
	   u(i,1,k)=0.d0 
	   v(i,1,k)=0.d0 
	   w(i,1,k)=0.d0
	   T(i,1,k)=T(i,2,k) 
	  enddo 
	  enddo 
	 endif 
!-------j=ny free-stream ----------
     if(npy .eq. npy0-1) then 
	 do k=1,nz 
	 do i=1,nx 
	  d(i,ny,k)=1.d0 
	  u(i,ny,k)=1.d0 
	  v(i,ny,k)=0.d0 
	  w(i,ny,k)=0.d0 
	  T(i,ny,k)=1.d0 
	 enddo 
	 enddo 
	endif 
!----------------------------------------
    end 	
