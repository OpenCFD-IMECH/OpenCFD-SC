
!=============================================================================================        
! Steger-Warming Splitting,   See Dexun Fu, Computational Aerodynamics, page 160-162 (in Chinese)
! Ref: Steger J I, warming R F, JCP 40, 263-293, 1981
  subroutine split_Stager_Warming(nx1,LAP1,c1,d1,u1,v1,w1,A1,A2,A3,fP,fm)
   Use flow_para
   implicit none
   integer nx1,LAP1,i
   real(kind=OCFD_REAL_KIND),dimension (1-LAP1:nx1+LAP1):: c1,d1,u1,v1,w1,A1,A2,A3
   real(kind=OCFD_REAL_KIND),dimension(1-LAP1:nx1+LAP1,5):: fp, fm 

    real(kind=OCFD_REAL_KIND):: ss,ss1,ak1,ak2,ak3,tmp0,tmp1,tmp2,tmp3,   &
          E1,E2,E3,E1P,E2P,E3P,E1M,E2M,E3M,   &
          vs,uc1,uc2,vc1,vc2,wc1,wc2,vvc1,vvc2,vv,W2,P2 
        
      tmp1=2.d0*(Para%gamma-1.d0)
      tmp2=1.d0/(2.d0*Para%gamma)
      tmp3=(3.d0-Para%gamma)/(2.d0*(Para%gamma-1.d0)) 

     do i=1-LAP1,nx1+LAP1
      ss=sqrt(A1(i)*A1(i)+A2(i)*A2(i)+A3(i)*A3(i))
      ss1=1.d0/ss
      ak1=A1(i)*ss1
      ak2=A2(i)*ss1
      ak3=A3(i)*ss1
      vs=A1(i)*u1(i)+A2(i)*v1(i)+A3(i)*w1(i)  
!c E1 is lamda1, lamda2 and lamda3; E2 is lamda4; E3 is lamda5 
      E1=vs
      E2=vs-c1(i)*ss
      E3=vs+c1(i)*ss
!----------------------------------------
      E1P=(E1+abs(E1))*0.5d0
      E2P=(E2+abs(E2))*0.5d0
      E3P=(E3+abs(E3))*0.5d0
 
      E1M=E1-E1P
      E2M=E2-E2P
      E3M=E3-E3P
!----------------------------------------
      tmp0=d1(i)/(2.d0*Para%gamma) 
      uc1=u1(i)-c1(i)*ak1
      uc2=u1(i)+c1(i)*ak1
      vc1=v1(i)-c1(i)*ak2
      vc2=v1(i)+c1(i)*ak2
      wc1=w1(i)-c1(i)*ak3
      wc2=w1(i)+c1(i)*ak3
      vvc1=(uc1*uc1+vc1*vc1+wc1*wc1)*0.5d0
      vvc2=(uc2*uc2+vc2*vc2+wc2*wc2)*0.5d0
      vv=(Para%gamma-1.d0)*(u1(i)*u1(i)+v1(i)*v1(i) +w1(i)*w1(i))   
      W2=tmp3*c1(i)*c1(i)
! --------The equation seems wrong ! P2 should be zero !
!      P2=tmp1*d1(i)*ak1*(ak2*w1(i)-ak3*v1(i))    ! ???????  Error
!--------------------------------------------------------
      fp(i,1)=tmp0*(tmp1*E1P+E2P+E3P)
      fp(i,2)=tmp0*(tmp1*E1P*u1(i)+E2P*uc1+E3P*uc2)
      fp(i,3)=tmp0*(tmp1*E1P*v1(i)+E2P*vc1+E3P*vc2)
      fp(i,4)=tmp0*(tmp1*E1P*w1(i)+E2P*wc1+E3P*wc2)
!     P2 should be zero !!!
!     fp(i,5)=tmp0*(E1P*vv+E2p*vvc1+E3P*vvc2+W2*(E2P+E3P)+P2*E1P)  ! Error 
      fp(i,5)=tmp0*(E1P*vv+E2p*vvc1+E3P*vvc2+W2*(E2P+E3P))
        
      fm(i,1)=tmp0*(tmp1*E1M+E2M+E3M)
      fm(i,2)=tmp0*(tmp1*E1M*u1(i)+E2M*uc1+E3M*uc2)
      fm(i,3)=tmp0*(tmp1*E1M*v1(i)+E2M*vc1+E3M*vc2)
      fm(i,4)=tmp0*(tmp1*E1M*w1(i)+E2M*wc1+E3M*wc2)

!     fm(i,5)=tmp0*(E1M*vv+E2M*vvc1+E3M*vvc2+W2*(E2M+E3M)+P2*E1M)  ! Error 
      fm(i,5)=tmp0*(E1M*vv+E2M*vvc1+E3M*vvc2+W2*(E2M+E3M))
     enddo
   end
!--------------------------------------------

 subroutine split_Local_LaxFriedrichs(nx1,LAP1,c1,d1,u1,v1,w1,A1,A2,A3,fP,fm)
   Use flow_para
   implicit none
   integer nx1,LAP1,i
   real(kind=OCFD_REAL_KIND),dimension (1-LAP1:nx1+LAP1):: c1,d1,u1,v1,w1,A1,A2,A3
   real(kind=OCFD_REAL_KIND),dimension(1-LAP1:nx1+LAP1,5):: fp, fm 

    real(kind=OCFD_REAL_KIND):: ss,un,p1,lamda0,E0,h0,tmp0,tmp1
     tmp0=1.d0/Para%gamma	
     tmp1=1.d0/(Para%gamma-1.d0)
	 
     do i=1-LAP1,nx1+LAP1
      ss=sqrt(A1(i)*A1(i)+A2(i)*A2(i)+A3(i)*A3(i))
	  un=A1(i)*u1(i)+A2(i)*v1(i)+A3(i)*w1(i)
	  lamda0=abs(un)+ss*c1(i)   
	  p1=tmp0*d1(i)*c1(i)*c1(i)
	  E0=tmp1*p1+0.5d0*d1(i)*(u1(i)*u1(i)+v1(i)*v1(i)+w1(i)*w1(i)) 
	  h0=E0+p1                     ! Total H (=E+p)
	  
	  fp(i,1)=0.5d0*(d1(i)*un +lamda0*d1(i))
	  fp(i,2)=0.5d0*(d1(i)*u1(i)*un+A1(i)*p1 + lamda0*d1(i)*u1(i))
	  fp(i,3)=0.5d0*(d1(i)*v1(i)*un+A2(i)*p1 + lamda0*d1(i)*v1(i))
	  fp(i,4)=0.5d0*(d1(i)*w1(i)*un+A3(i)*p1 + lamda0*d1(i)*w1(i))
	  fp(i,5)=0.5d0*(h0*un + lamda0*E0)	  

	  fm(i,1)=0.5d0*(d1(i)*un -lamda0*d1(i))
	  fm(i,2)=0.5d0*(d1(i)*u1(i)*un+A1(i)*p1 - lamda0*d1(i)*u1(i))
	  fm(i,3)=0.5d0*(d1(i)*v1(i)*un+A2(i)*p1 - lamda0*d1(i)*v1(i))
	  fm(i,4)=0.5d0*(d1(i)*w1(i)*un+A3(i)*p1 - lamda0*d1(i)*w1(i))
	  fm(i,5)=0.5d0*(h0*un - lamda0*E0)	  	  
	 enddo 
	 end 