
! interpolation for compression-corner type mesh  ( mesh uniform in Z- direction )

    implicit none   
    integer,parameter:: Periodic=1,Non_Periodic=0
    integer:: Nvar,nx1,ny1,nz1,nx2,ny2,nz2
	real*8,allocatable,dimension(:):: Z1, Z2, W1, W2  
    real*8,allocatable,dimension(:,:) :: x1,y1,x2,y2,Ak
	real*8,allocatable,dimension(:,:,:) :: V1, V2, U2
    real*8,allocatable,dimension(:,:,:,:):: U1,Un
    real*8:: ZL,tt   
    integer,allocatable,dimension(:) ::  ks_k	
    integer:: Istep ,i,j,k,m
	
	print*, "3D Interploation  for  compression corner ..."
	print*, "please input nx1,ny1,nz1  (orginal grid)"
	read(*,*) nx1, ny1, nz1 
	print*, "please input nx2, ny2, nz2 (new grid)"
	read(*,*) nx2, ny2, nz2 
	Nvar=5
	ZL=1.d0 
	
   allocate( z1(nx1), z2(nz2),  ks_k(nz2) , Ak(6,nz2), W1(nz1), W2(nz2))
   allocate(x1(nx1,ny1), y1(nx1,ny1), x2(nx2,ny2), y2(nx2,ny2))
   allocate(V1(nx1,ny1,Nvar),V2(nx2,ny2,Nvar))
   allocate(U1(nx1,ny1,nz1,Nvar),Un(nx2,ny2,nz1,Nvar),U2(nx2,ny2,nz2) )
   
!--------read coordinate ------------------------------   
   open(99,file="mesh-old.dat",form="unformatted") 
   read(99) x1
   read(99) y1 
   close(99) 
   
   open(100,file="mesh-new.dat",form="unformatted")
   read(100) x2 
   read(100) y2 
   close(100) 

!-------read flow data ----------------
   print*, "read old data ..."
   open(99,file="data-old.dat",form="unformatted")
   read(99) Istep, tt 
   print*, "Istep, tt=", Istep, tt
   do m=1,Nvar
    print*, "read 3d data ..." 
	do k=1,nz1
	read(99) U1(:,:,k,m) 
	enddo 
   enddo 
   close(99) 
  
!------------Inter ploation in x-y plane --------------------  
  print*, "read flow data OK , Interploation in x-y plane ..."
  
  do k=1,nz1    
   print*, "k=", k
    do m=1,Nvar
	do j=1,ny1
	do i=1,nx1 
	 V1(i,j,m)=U1(i,j,k,m)
	enddo
	enddo
	enddo 
    
	call interploation_curve2d(Nvar,nx1,ny1,nx2,ny2,x1,y1,V1,x2,y2,V2)
  
   do m=1,Nvar 
   do j=1,ny2 
   do i=1,nx2 
    Un(i,j,k,m)=V2(i,j,m) 
   enddo
   enddo
   enddo 
  enddo
 
  open(99,file="flow-z1-old.dat")
  write(99,*) "variables=x,y,T"
  write(99,*) "zone i=" ,nx1, "j=", ny1 
  do j=1,ny1 
  do i=1,nx1
  write(99,"(3F16.8)") x1(i,j), y1(i,j), U1(i,j,1,5)
  enddo
  enddo 
  
  open(99,file="flow-z1-new.dat")
  write(99,*) "variables=x,y,T"
  write(99,*) "zone i=" ,nx2, "j=", ny2 
  do j=1,ny2 
  do i=1,nx2
  write(99,"(3F16.8)") x2(i,j), y2(i,j), Un(i,j,1,5)
  enddo
  enddo 
  
!-----------------------------------------------------------
!--interploation in z- direction ----------------

  	  do k=1,nz1
	   z1(k)=(k-1.d0)/nz1
	  enddo
	  do k=1,nz2
	   z2(k)=(k-1.d0)/nz2 
	  enddo
      call get_ks_Ai(nz1,nz2,z1,z2,ks_k,Ak,Periodic, ZL) 

      open(55,file="data-new.dat",form="unformatted") 
	  write(55) Istep, tt

	  do m=1,Nvar
        print*, "interpolation-Z... ",m
 
!---------------interpolation in z- direction
         do j=1,ny2
         do i=1,nx2
          w1=Un(i,j,1:nz1,m)           
	      call interpolation6(nz1,nz2,ks_k,Ak,w1,w2,Periodic)
	      U2(i,j,1:nz2)=w2
         enddo
		 enddo

!------------------------------------------          
 	   call write3d(55,U2,nx2,ny2,nz2)
 	   enddo
      close(55)	  

!----------test in y-z-direction	 
	  i=nx2/2
      open(44,file="flow2d-yz-1.dat")
      write(44,*) "variables=x,y,T"
      write(44,*) "zone i= ",ny2 ,  " j= ", nz1
	  do k=1,nz1
	  do j=1,ny2
	   write(44,"(3f15.6)") (j-1.d0)/ny2, (k-1.0)/nz1, Un(i,j,k,5)
	  enddo
	  enddo
      close(44)

       open(44,file="flow2d-yz-2.dat")
        write(44,*) "variables=x,y,T"
        write(44,*) "zone i= ", ny2,  " j= ", nz2
	    do k=1,nz2
	    do j=1,ny2
	     write(44,"(3f15.6)") (j-1.d0)/ny2, (k-1.0)/nz2,U2(i,j,k)
	    enddo
	    enddo
        close(44)
!=======================================

	deallocate( Z1, Z2, W1, W2,x1,y1,x2,y2,Ak, V1, V2, U1,U2,Un, Ks_k)

   end


!-----------2nd order interploation ----------------------
! Nvar : Number of variables
 subroutine interploation_curve2d(Nvar,nx1,ny1,nx2,ny2,x1,y1,V1,x2,y2,V2)
 implicit none
  integer::  Nvar,nx1,ny1,nx2,ny2
  integer::  i,j,i0,j0,m
  real*8::   x1(nx1,ny1),y1(nx1,ny1),x2(nx2,ny2),y2(nx2,ny2),V1(nx1,ny1,Nvar),V2(nx2,ny2,Nvar) 
  real*8,allocatable,dimension(:,:,:):: Vx,Vy
  integer,allocatable,dimension(:,:) :: ist, jst        ! the nearest (i,j)
	print*, "in inter_ploation 2D ..."

!----------------------------
  allocate(ist(nx2,ny2), jst(nx2,ny2))
  allocate(Vx(nx1,ny1,Nvar),Vy(nx1,ny1,Nvar)) 
 
 call  comput_grident(Nvar,nx1,ny1,x1,y1,V1,Vx,Vy)

 call find_nearest_ij(nx1,ny1, nx2,ny2, x1,y1,x2,y2, ist, jst)                  ! Find the nearest point

!-------------Taylor expansion ----------------------------------------
 do j=1,  ny2
 do i=1,  nx2
    i0=ist(i,j)                  ! 最近的点坐标
	j0=jst(i,j) 
 do m=1,Nvar 
   V2(i,j,m)= V1(i0,j0,m) + Vx(i0,j0,m)*(x2(i,j)-x1(i0,j0))+ Vy(i0,j0,m)*(y2(i,j)-y1(i0,j0))
 enddo 
      
 enddo
 enddo
 
 deallocate(ist, jst, Vx, vy)
  end


!------------计算梯度---------------------------------------------
   subroutine  comput_grident(Nvar,nx,ny, xx, yy, V1, Vx,Vy)
   implicit none
   integer:: Nvar,nx,ny,i,j
   real*8:: xx(nx,ny), yy(nx,ny),V1(nx,ny,Nvar),Vx(nx,ny,Nvar),Vy(nx,ny,Nvar)
   real*8::  xi, xj, yi, yj,  Aix, Aiy, Ajx, Ajy, Ajac,  Vi(Nvar),Vj(Nvar)

   do j=1,ny
   do i=1,nx
    if(i .eq. 1) then 
	 xi=xx(2,j)-xx(1,j)
	 yi=yy(2,j)-yy(1,j) 
	 Vi(:)=V1(2,j,:)-V1(1,j,:)
	else if(i .eq. nx )then 
     xi=xx(nx,j)-xx(nx-1,j) 
     yi=yy(nx,j)-yy(nx-1,j) 
     Vi(:)=	V1(nx,j,:)-V1(nx-1,j,:) 
	else 
     xi=(xx(i+1,j)-xx(i-1,j))*0.5d0 
     yi=(yy(i+1,j)-yy(i-1,j))*0.5d0 
     Vi(:)=	(V1(i+1,j,:)-V1(i-1,j,:))*0.5d0  	 
	endif 

    if(j .eq. 1) then 
	 xj=xx(i,2)-xx(i,1)
	 yj=yy(i,2)-yy(i,1) 
	 Vj(:)=V1(i,2,:)-V1(i,1,:)
	else if(j .eq. ny )then 
     xj=xx(i,ny)-xx(i,ny-1) 
     yj=yy(i,ny)-yy(i,ny-1) 
     Vj(:)=	V1(i,ny,:)-V1(i,ny-1,:) 
	else 
     xj=(xx(i,j+1)-xx(i,j-1))*0.5d0 
     yj=(yy(i,j+1)-yy(i,j-1))*0.5d0 
     Vj(:)=	(V1(i,j+1,:)-V1(i,j-1,:))*0.5d0  	 
	endif 	
	
	  Ajac=1.d0/(xi*yj-yi*xj)
	  Aix=yj*Ajac ;  Aiy=-xj*Ajac;  Ajx=-yi*Ajac ; Ajy=xi*Ajac                     ! Jacobian
	  Vx(i,j,:)=Vi(:)*Aix+Vj(:)*Ajx 
	  Vy(i,j,:)=Vi(:)*Aiy+Vj(:)*Ajy

  enddo
  enddo
  


  end

!-----------------------搜索--------------------------------------------------------------------------------
! (xn,yn) old coordinate;   (x2,y2) new coordinate
   subroutine find_nearest_ij(nx0,ny0, nx2,ny2, xn, yn, x2, y2, ist, jst)
   implicit none
   integer:: nx0,ny0, nx2, ny2, i, j, i1,j1, ia,ib, ja, jb, Iflag_Half, i0, j0
   real*8:: xn(nx0,ny0), yn(nx0,ny0), x2(nx2,ny2), y2(nx2,ny2)
   integer::  ist(nx2,ny2), jst(nx2,ny2)
   real*8:: d0, dd
   integer, parameter:: LP=10
 !--------------------------------------------------
   do  j1=1,  ny2
   do  i1=1,  nx2
    if(i1 .ne. 1)  then
	   ia=max(ist(i1-1, j1)-LP, 1) ;    ib=min(ist(i1-1, j1)+LP, nx0)
	   ja=max(jst(i1-1, j1)-LP, 1) ;    jb=min(jst(i1-1, j1)+LP, ny0)
    else
	  if(j1 .ne. 1) then
	   ia=max(ist(i1, j1-1)-LP, 1) ;    ib=min(ist(i1, j1-1)+LP, nx0)
	   ja=max(jst(i1, j1-1)-LP, 1) ;    jb=min(jst(i1, j1-1)+LP, ny0)
	 else
       ia=1;  ib=nx0
	   ja=1;  jb=ny0
	 endif
   endif
       
   i0=ia;  j0=ja ;  d0=1.d20             ! d0, a larage number
   do j=ja, jb
   do i=ia, ib
      dd=(x2(i1,j1)-xn(i,j))**2+(y2(i1,j1)-yn(i,j))**2
	  
	 if(dd < d0) then
	    d0=dd
		i0=i;   j0=j
	 endif
	enddo
	enddo
	ist(i1,j1)=i0
	jst(i1,j1)=j0

  enddo
  enddo
  end	 	    
  
  
!-------------------------------------------------------
     subroutine get_ks_Ai(nx1,nx2,xx1,xx2,ks_i,Ai,IF_Periodic, ZL_Periodic)
	 implicit doubleprecision(a-h,o-z)
	 integer,parameter:: LAP=6
	 integer nx1,nx2,ks_i(nx2)
	 real*8 xx1(nx1),xx2(nx2),Ai(6,nx2), ZL_Periodic
	 real*8:: xx0(1-LAP:nx1+LAP)    ! periodical spaned x

	 do j=1,nx1
	   xx0(j)=xx1(j)
	 enddo
	 do j=1-LAP, 0
	  xx0(j)=xx1(nx1+j)-ZL_Periodic
	enddo
	do j=nx1+1, nx1+LAP
	 xx0(j)=xx1(j-nx1)+ZL_Periodic
   enddo
    print*, "xx0=", xx0

	 do j=1,nx2
        if(xx2(j) .lt. xx1(1) ) then
         ks_i(j)=0
         goto 100
	    endif

	  do i=1,nx1-1
	   if(xx2(j) .ge. xx1(i) .and. xx2(j) .lt. xx1(i+1)) then
	     ks_i(j)=i
	     goto 100
	   endif
	  enddo
         ks_i(j)=nx1
100     continue
       enddo

	    Ai=0.d0


         do i=1,nx2
		     ka=4-ks_i(i)
		  if(ka .lt. 1 .or. IF_Periodic .eq. 1) ka=1
             kb=nx1+3-ks_i(i)
          if(kb .gt. 6 .or. IF_Periodic .eq. 1) kb=6

		do k=ka,kb
!	       ik=MOD_N(ks_i(i)+k-3,nx1)
	 	   ik=ks_i(i)+k-3
          Ai(k,i)=1.d0
	   do km=ka,kb
	      ikm=ks_i(i)+km-3
!	     ikm=MOD_N(ks_i(i)+km-3,nx1)

	   if(km .ne. k) then
!	   Ai(k,i)=Ai(k,i)*(xx2(i)-xx1(ikm))/(xx1(ik)-xx1(ikm))
	   Ai(k,i)=Ai(k,i)*(xx2(i)-xx0(ikm))/(xx0(ik)-xx0(ikm))
	   endif
	   enddo
	   enddo

     enddo


       end

!----------------------------------------

	  subroutine interpolation6(nx1,nx2,ks_i,Ai,ux1,ux2,IF_Periodic)
	  implicit doubleprecision(a-h,o-z)
	  integer nx1,nx2,ks_i(nx2)
	  real*8 Ai(6,nx2),ux1(nx1),ux2(nx2)

      do i=1,nx2
	  
	    ka=4-ks_i(i)
	    if(ka .lt. 1 .or. IF_Periodic .eq. 1 ) ka=1
        kb=nx1+3-ks_i(i)
        if(kb .gt. 6  .or. IF_Periodic .eq. 1) kb=6

        ux2(i)=0.d0
        do k=ka,kb
! 	     ik=ks_i(i)+k-3
 	     ik=MOD_N(ks_i(i)+k-3,nx1)
         ux2(i)=ux2(i)+Ai(k,i)*ux1(ik)
        enddo

         if(ks_i(i) .eq. 0 .and. IF_Periodic .eq. 0) ux2(i)=ux1(1)     ! do not ext-interpolate
         if(ks_i(i) .eq. nx1 .and. IF_Periodic .eq. 0) ux2(i)=ux1(nx1)

       enddo
       end

!---------------------------------------
    function MOD_N(k,N)
	 integer k,N,MOD_N
	 MOD_N=k
	 if(k .lt. 1) MOD_N=k+N
	 if(K .gt. N) MOD_N=k-N
	end

!----------------------------------------------
     subroutine read3d(no,U,nx,ny,nz)
       implicit doubleprecision (a-h,o-z)
       real*8 U(nx,ny,nz)
       do k=1,nz
       read(no) U(:,:,k)
       enddo
     end
!----------------------------------------------
     subroutine write3d(no,U,nx,ny,nz)
       implicit doubleprecision (a-h,o-z)
       real*8 U(nx,ny,nz)
       do k=1,nz
       write(no) U(:,:,k)
       enddo
     end  