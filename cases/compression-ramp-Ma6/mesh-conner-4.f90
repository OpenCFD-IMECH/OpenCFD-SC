! Mesh generater of a compression ramp V3.3
! Copy right by Lixinliang 
! Ver 3.0: 2008-8 
! Ver 3.1: 添加了边界正交化模块
! ver 3.2: allocatable array used in Jocabian() 
! Ver 3.3: 2020-4-26;  到壁面距离为变数；  （x<0 时维持在hy_wall1,  x>0 时线性变化到hy_wall2
! Ver 4:  2021-3-12: data for OpenCFD 2 (OCFD-grid.dat)
     implicit doubleprecision (a-h,o-z)
     real*8,parameter:: PI=3.1415926535897932     
     real*8,allocatable:: xa(:),xb(:), ya(:),yb(:),sx(:),sy(:),xx(:,:),yy(:,:),zz(:),hy_wall(:)
	 real*8 XL,Lz
	 open(66,file='grid-conner.in')
	 read(66,*)
	 read(66,*)
	 read(66,*) nx,ny,nz
     read(66,*)
	 read(66,*) alfa,High,hy_wall1,hy_wall2,Lz	 
	 read(66,*)
     read(66,*) XL_conner,nx_conner,XL_inlet,nx_inlet,dx_inlet,nx_buff,alfax_buff
	 close(66)

 	 allocate(xa(nx),xb(nx) ,ya(nx),yb(nx),sx(nx),sy(ny),xx(nx,ny),yy(nx,ny),hy_wall(nx),zz(nz))
 
!c-------------------------------
      if(nx_conner+nx_inlet+nx_buff .ne. nx) then
	    print*, "Error ! nx_conner+nx_inlet+nx_buff != nx" 
		stop
	  endif

	  alfa=alfa*PI/180.d0
	  hx0=2.d0*XL_conner/(nx_conner-1.d0)
      R=XL_conner/tan(alfa/2.d0)-High

  call   grid1d_transition(nx_inlet+1,xa, XL_inlet,dx_inlet, hx0)
  	  do i=1,nx_inlet
      xa(i)=xa(i)-XL_conner-XL_inlet
	  ya(i)=0.d0
	  xb(i)=xa(i)
	  yb(i)=High
      enddo

!------conner region ----------------------
      i0=nx_inlet
      do i=1,nx_conner
	  seta=alfa*(i-1.d0)/(nx_conner-1.d0)
	  if(seta .le. alfa/2.d0) then
	   xa(i+i0)=  (i-1.d0)*hx0-XL_conner
	   ya(i+i0)=  0.d0
	  else
	   xa(i+i0)=  ((i-1.d0)*hx0-XL_conner)*cos(alfa)
	   ya(i+i0)=  ((i-1.d0)*hx0-XL_conner)*sin(alfa)
      endif
       xb(i+i0)=R*sin(seta)-XL_conner
	   yb(i+i0)=R*(1.d0-cos(seta))+High
	 enddo
!----buff region -------------------------
    i0=nx_inlet+nx_conner
	do i=1,nx_buff 
	dx=(xa(i+i0-1)-xa(i+i0-2))*alfax_buff
    xa(i+i0)=xa(i+i0-1)+dx
	ya(i+i0)=ya(i+i0-1)+dx*tan(alfa)
	xb(i+i0)=xb(i+i0-1)+dx
	yb(i+i0)=yb(i+i0-1)+dx*tan(alfa)
	enddo
!-------------------------------------------

     open(33,file='grid1d.dat')
     do i=1,nx
	 write(33,'(i6, 6f16.8)') i, xa(i),ya(i),xb(i),yb(i)
	 enddo
!------Computing hy in wall -----------(hy_wall=hy_wall1 in x<0 regin, keep linear to hy_wall2 in x>0 regin)
      do i=1,nx
	  if(xa(i)<=0) then 
	   hy_wall(i)=hy_wall1
	  else
	   hy_wall(i)=hy_wall1+xa(i)/xa(nx)*(hy_wall2-hy_wall1) 
	  endif
	 enddo
!----------------------------------------------

	    
	 do i=1,nx
	  SL=sqrt((xb(i)-xa(i))**2+(yb(i)-ya(i))**2)
	  call getsy(ny,sy,SL,hy_wall(i))
	  do j=1,ny
	  xx(i,j)=xa(i)+(xb(i)-xa(i))*sy(j)
	  yy(i,j)=ya(i)+(yb(i)-ya(i))*sy(j)
	  enddo
	 enddo
        call Wall_Orthonormalization(nx,ny,xx,yy)

        open(44,file="y1d.dat")
        do j=1,ny
        write(44,*) j, yy(1,j)
        enddo
        close(44) 
		
		do k=1,nz 
		 zz(k)=Lz*(k-1.d0)/nz
		enddo 
		open(99,file="OCFD-grid.dat",form="unformatted")
		write(99) xx 
		write(99) yy 
		write(99) zz 
		close(99) 
       open(33,file='grid.dat')
	   write(33,*) 'variables=x,y'
	   write(33,*) 'zone i=',nx , 'j=',ny
	   do j=1,ny
	   do i=1,nx
	   write(33,'(2f15.6)') xx(i,j),yy(i,j)
	   enddo
	   enddo
		
		
	 deallocate(xa,xb ,ya,yb,sx,sy,xx,yy,zz,hy_wall)	
	 end

!c=================================================================
         subroutine getsy(ny,sy,SL,deltx)
 	     implicit doubleprecision (a-h,o-z)
	     real*8 sy(ny)
         real*8,save:: b=3.5d0
	     dy=1.d0/(ny-1.d0)
!---------------------------------------
         delta=deltx/SL
 ! using Newton method get coefficient
 100     continue
!        fb=exp(b/(ny-1.d0))-1.d0-delta*(exp(b)-1.d0)
!	     fbx=exp(b/(ny-1.d0))/(ny-1.d0)-delta*exp(b)
         fb=(exp(b/(ny-1.d0))-1.d0)/(exp(b)-1.d0)-delta
         fbx=(exp(b/(ny-1.d0))/(ny-1.d0)*(exp(b)-1.d0) -      &     
         (exp(b/(ny-1.d0))-1.d0)*exp(b)  )/((exp(b)-1.d0))**2
         bnew=b-fb/fbx
         if(abs(b-bnew) .gt. 1.d-6) then
	     b=bnew
!           print*, "b=",b
		 goto 100
	     endif
	    
         b=bnew 
         a=1.d0/(exp(b)-1.d0)
         do j=1,ny
         s=(j-1.d0)*dy
         sy(j)=a*(exp(s*b)-1.d0)
         enddo
       end




!c======== x direction ------------------
       subroutine getsx(nx,nx_buff,sx,alfax,alfax_buff)
       implicit doubleprecision (a-h,o-z)
	   dimension sx(nx)
         nx1=nx-nx_buff
	   dx=1.d0/(nx1-1.d0)
	   A=1.d0/(exp(alfax)-1.d0)
         do i=1,nx1
         s=(i-1./2.)*dx
	   sx(i)=A*(exp(s*alfax)-1.d0)
         enddo
	   do i=nx1+1,nx
	   sx(i)=sx(i-1)+alfax_buff*(sx(i-1)-sx(i-2))
	   enddo
         open(55,file="sx.dat")
	   do k=1,nx-1
	   write(55,*) k, sx(k),sx(k+1)-sx(k)
	   enddo
	  end





!--------------------------------------------------------------
! 模块2: 过渡网格，采用三次函数过渡，光滑连接两端
! 接口：nx 网格数； Length 长度； dx1 首端网格间距； dx2 尾端网格间距
! 通常要求 平均网格间距 Length/(nx-1.) 介于 dx1 与dx2之间
 
    subroutine grid1d_transition(nx,xx, Length, dx1,dx2)
    implicit doubleprecision (a-h,o-z)
	real*8 Length, xx(nx)
      dx0= Length/(nx-1.d0)  
      if( dx0 .gt. dmax1(dx1,dx2) .or. dx0 .lt. dmin1(dx1,dx2) ) then
        print*, "warning !  dx0 should between dx1 and dx2 !!!"
        print*, "dx0  dx1 and dx2 =", dx0, dx1, dx2
      endif

    do i=1,nx
      x=(i-1.d0)/(nx-1.d0)
!     Ah0=(1.d0+2.d0*x)*(x-1.d0)**2    ! 三次Hermit插值的4个基函数
      Ah1=(1.d0-2.d0*(x-1.d0))*x*x
      Sh0=x*(x-1.d0)**2
      Sh1=(x-1.d0)*x*x
      xx(i)=Length*Ah1+dx1*(nx-1.d0)*Sh0+dx2*(nx-1.d0)*Sh1
    enddo
   end

!----------------------------------------------------------
! See "计算空气动力学 p301 or Obayashi et al AIAA-86-0338,1986
     subroutine Wall_Orthonormalization(nx,ny,xx,yy)
     implicit doubleprecision (a-h,o-z)
	 dimension xx(nx,ny),yy(nx,ny)
	 real*8 dx(nx),dy(nx),dx1(nx),dy1(nx)

     parameter(Maxstep=10)
!	 parameter (b=3.d0) 
     parameter (b=2.d0)          ! 2020-4-26
     Jmax=20

	 do mstep= 1, Maxstep
	 do j=2,Jmax
     dx(1)=0.d0; dy(1)=0.d0; dx(nx)=0.d0; dy(nx)=0.d0;
	 do i=2,nx-1
	 qn=sqrt((xx(i+1,j)-xx(i-1,j))**2+(yy(i+1,j)-yy(i-1,j))**2)
	 qnx=-(yy(i+1,j)-yy(i-1,j))/qn
	 qny= (xx(i+1,j)-xx(i-1,j))/qn
     s= (xx(i,j)-xx(i,j-1))*qnx+(yy(i,j)-yy(i,j-1))*qny
	 dx(i)=s*qnx-(xx(i,j)-xx(i,j-1))
	 dy(i)=s*qny-(yy(i,j)-yy(i,j-1))
     enddo

	 do i=2,nx-1
	  dx1(i)=(dx(i-1)+2.d0*dx(i)+dx(i+1))/4.d0
	  dy1(i)=(dy(i-1)+2.d0*dy(i)+dy(i+1))/4.d0
	  xx(i,j)=xx(i,j)+exp(-b*j/Jmax)/2.d0*dx1(i)
	  yy(i,j)=yy(i,j)+exp(-b*j/Jmax)/2.d0*dy1(i)
	 enddo
     enddo
	enddo

	end



