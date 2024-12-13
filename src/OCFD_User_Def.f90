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
! �û����ж��壨�������ĸ�ʽ
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
! �û���ʽ��������ܵ㣬  ������Զ���ĸ�ʽ�ṩ
! [Ka_P,Kb_P] Ϊ��ͨ����ʽ��ʹ�õ���������� 
! ����WENOP ��ʽ�� ����h(i+1/2)ʱ�� ʹ���� i-2, i-1, i, i+1, i+2 �㣬 ���Ka_P=-2, Kb_P=2 

! [Ka_M,Kb_M] Ϊ��ͨ����ʽ��ʹ�õ����������	
! ����WENOM ��ʽ�� ����h(i+1/2)ʱ�� ʹ���� i-1, i, i+1, i+2, i+3 �㣬 ���Ka_M=-1, Kb_M=3 
 
	   subroutine Stencil_Scheme_User(Ka_P, Kb_P,  Ka_M, Kb_M)
	   implicit none 
       integer:: Ka_P, Kb_P, Ka_M, Kb_M 	   
	   

!------------������WENO ��ʽ�����ӣ� �밴���Լ��ĸ�ʽ�޸�   		
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
! �û��Զ���ĸ�ʽ��  ��ֵͨ�� ����ͨ����

! Ka, Kb ��������ܵ�������   ����i+1/2��ֵͨ��ʱ�� ʹ�õ����������Ϊ[i+Ka, i+Kb]

! ����WENOP ��ʽ�� ����h(i+1/2)ʱ�� ʹ���� i-2, i-1, i, i+1, i+2 �㣬 ���Ka=-2, Kb=2 

	subroutine hh_Scheme_USER_P(Ka,Kb,v,hh)           
     Use OCFD_constants
     implicit none
     integer Ka,Kb
     real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
!      Add you scheme code here !  ����Լ��Ĵ��� �������� WENO5P ��ʽ�����ӣ�
 
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

!  ��ֵͨ�� ����ͨ����   
    subroutine hh_Scheme_USER_M(Ka,Kb,v,hh)            
     Use OCFD_constants
     implicit none
     integer Ka,Kb
     real(kind=OCFD_REAL_KIND):: v(Ka:Kb), hh
 
!      Add you scheme code here !  ����Լ��Ĵ��� �������� WENO5M ��ʽ�����ӣ�
 
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











