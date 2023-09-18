!-------------User defined code ---------------------------------
! ocfd_bc_User()  :   boundary conditions defined by user
! ana_user():         post analysis code defined by user 
! hh_Scheme_USER_P() :  scheme (positive flux) defined by user 
! hh_Scheme_USER_M() :  Scheme (Negative flux) defined by user 
! Stencil_Scheme_User(): use defined Scheme's Stencil 




! boundary conditions defined by user---------------------
	subroutine ocfd_bc_User
    use flow_data
	implicit none 
!--------Input your boundary condition code here !	
	print*, "Input your boundary condition code here !"
!      input your boundary condition code here !
    end 	
	
	
!-------------------------------------------------------------
!  Post analysis code defined by user     
   subroutine ana_user(ana_data)  
   Use flow_data     
   implicit none
   real(kind=OCFD_REAL_KIND):: ana_data(50)  
!      input your post-analysis code here !!!   
   end  



! User's Scheme 
! 用户自行定义（开发）的格式
!------------------------------------------------------     
	 subroutine OCFD_Scheme_USER_P(v,hh,nx,LAP,ib,ie)
     Use OCFD_constants
	 Use Scheme_para
     implicit none
     integer nx,LAP,i,Ka,Kb,ib,ie
     real(kind=OCFD_REAL_KIND):: v(1-LAP:nx+LAP), hh(0:nx)
     Ka=Scheme%Ka1 
	 Kb=Scheme%Kb1
!	 do i=0,nx 
     do i=ib,ie
	  call hh_Scheme_USER_P(Ka,Kb,v(i+Ka:i+Kb),hh(i))
	 enddo 
	 end 
	 
	 subroutine OCFD_Scheme_USER_M(v,hh,nx,LAP,ib,ie)
     Use OCFD_constants
	 Use Scheme_para
     implicit none
     integer nx,LAP,i,Ka,Kb,ib,ie 
     real(kind=OCFD_REAL_KIND):: v(1-LAP:nx+LAP), hh(0:nx)
     Ka=Scheme%Ka2 
	 Kb=Scheme%Kb2   
!	 do i=0,nx 
     do i=ib,ie
	  call hh_Scheme_USER_M(Ka,Kb,v(i+Ka:i+Kb),hh(i))   
	 enddo 
	 end 	 
!----------------------------------------------------------------
! 用户格式的网格基架点，  请根据自定义的格式提供
! [Ka_P,Kb_P] 为正通量格式所使用的网格点区间 
! 例如WENOP 格式， 计算h(i+1/2)时， 使用了 i-2, i-1, i, i+1, i+2 点， 因此Ka_P=-2, Kb_P=2 

! [Ka_M,Kb_M] 为负通量格式所使用的网格点区间	
! 例如WENOM 格式， 计算h(i+1/2)时， 使用了 i-1, i, i+1, i+2, i+3 点， 因此Ka_M=-1, Kb_M=3 
 
	   subroutine Stencil_Scheme_User(Ka_P, Kb_P,  Ka_M, Kb_M)
	   implicit none 
       integer:: Ka_P, Kb_P, Ka_M, Kb_M 	   
	   

!------------下面是WENO 格式的例子， 请按照自己的格式修改   		
!    Add you define here !!! 

!       Ka_P=-2
!	    Kb_P=2
!	    Ka_M=-1
!	    Kb_M=3

        print*, "Please set Ka_P,Kb_P;  Ka_M, Kb_M for your scheme !"        
!        stop	   
!-----------------------------------------	   
      end 
		 
		   
		   
!-------Numerical flux for Positive flux   
! 用户自定义的格式：  数值通量 （正通量）

! Ka, Kb 是网格基架点描述；   计算i+1/2数值通量时， 使用的网格点区间为[i+Ka, i+Kb]

! 例如WENOP 格式， 计算h(i+1/2)时， 使用了 i-2, i-1, i, i+1, i+2 点， 因此Ka=-2, Kb=2 

	subroutine hh_Scheme_USER_P(Ka,Kb,v,hh)           
     Use OCFD_constants
     implicit none
     integer Ka,Kb
     real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
!      Add you scheme code here !  添加自己的代码 （下面是 WENO5P 格式的例子）
 
!    real(kind=OCFD_REAL_KIND)::  S0,S1,S2,a0,a1,a2,am,q03,q13,q23
!     real(kind=OCFD_REAL_KIND),parameter::   ep=1.d-6 ,   C03=3.d0/10.d0 ,    C13=3.d0/5.d0 ,    C23=1.d0/10.d0

!         S0=13.d0/12.d0*(v(0)-2.d0*v(1)+v(2))**2+  1.d0/4.d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
!         S1=13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+  1.d0/4.d0*(v(-1)-v(1))**2
!         S2=13.d0/12.d0*(v(-2)-2.d0*v(-1)+v(0))**2+  1.d0/4.d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
!         a0=C03/((ep+S0)**2)
!         a1=C13/((ep+S1)**2)
!         a2=C23/((ep+S2)**2)
!         am=a0+a1+a2
!         q03=1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
!         q13=-1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)
!         q23=1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
!         hh=(a0*q03+a1*q13+a2*q23)/am

   end                   

!  数值通量 （负通量）   
    subroutine hh_Scheme_USER_M(Ka,Kb,v,hh)            
     Use OCFD_constants
     implicit none
     integer Ka,Kb
     real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
 
!      Add you scheme code here !  添加自己的代码 （下面是 WENO5M 格式的例子）
 
!     real(kind=OCFD_REAL_KIND)::  S0,S1,S2,a0,a1,a2,am,q03,q13,q23
!     real(kind=OCFD_REAL_KIND),parameter::   ep=1.d-6 ,   C03=3.d0/10.d0 ,    C13=3.d0/5.d0 ,    C23=1.d0/10.d0  
 
!       S0=13.d0/12.d0*(v(1)-2.d0*v(0)+v(-1))**2+  1.d0/4.d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
!       S1=13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+  1.d0/4.d0*(v(2)-v(0))**2
!       S2=13.d0/12.d0*(v(3)-2.d0*v(2)+v(1))**2+  1.d0/4.d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
!       a0=C03/((ep+S0)**2)
!       a1=C13/((ep+S1)**2)
!       a2=C23/((ep+S2)**2)
!       am=a0+a1+a2
!       q03=1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)
!       q13=-1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)
!       q23=1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)
!       hh=(a0*q03+a1*q13+a2*q23)/am

    end











