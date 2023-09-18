! Generate initial data for isotropic turbulence (d,u,v,w,T)
! d=T=1;  u,v,w with energy spectrum  and random phase
       implicit none
        integer,parameter:: nx=128, ny=nx, nz=nx
        real*8, allocatable,dimension(:,:,:):: ur,ui,vr,vi,wr,wi
        integer:: i,j,k,k0,k1,k2,k3,kk
        real*8:: A, Ek0,tao,a1,a2,a3,v1r,v1i,v2r,v2i,v3r,v3i,vkr,vki,Er,Ei,sr,Ak
        real*8:: PAI=3.1415926535897932d0
!----------------------------------------------------        
        allocate(ur(nx,ny,nz),ui(nx,ny,nz),vr(nx,ny,nz),vi(nx,ny,nz),wr(nx,ny,nz),wi(nx,ny,nz))
        
 !--------------------------------------------------
 
           k0=8
          A=0.00013d0
          Ek0=3.d0*A/64.d0*sqrt(2.d0*PAI)*k0**5              ! Total energy
          tao=sqrt(32.d0/A)*(2.d0*PAI)**0.25d0* dble(k0)**(-7.d0/2.d0)   ! eddy turn-over time
          print*, 'tao=',tao,'EK0=',Ek0
!-------------------------------------------------
    
         do  k=1,nz
         do  j=1,ny
         do  i=1,nx
          k1=i-1
          k2=j-1
          k3=k-1
         if(k1 .ge. nx/2) k1=i-1-nx
         if(k2 .ge. ny/2) k2=j-1-ny
         if(k3 .ge. nz/2) k3=k-1-nz
      
       kk=(k1*k1+k2*k2+k3*k3)
       if(kk.ne.0) then
! ------- Energy spectrum --------------------
!  you can modify it      
        Ak=sqrt(2.d0/3.d0)*sqrt(A*kk*exp(-2.d0*kk/(k0*k0))/(4.d0*PAI))
!--------------------------------------------------
!       random phase
         call random_number(a1)
         call random_number(a2)
         call random_number(a3)
               
         a1=a1*2.d0*PAI
         a2=a2*2.d0*PAI
         a3=a3*2.d0*PAI
       
       v1r=Ak*sin(a1) ;  v1i=Ak*cos(a1)
       v2r=Ak*sin(a2) ;  v2i=Ak*cos(a2)
       v3r=Ak*sin(a3) ;  v3i=Ak*cos(a3)
       
	   vkr=k1*v1r+k2*v2r+k3*v3r
	   vki=k1*v1i+k2*v2i+k3*v3i
	 
	    ur(i,j,k)=k1*vkr/kk-v1r
	    ui(i,j,k)=k1*vki/kk-v1i
	    vr(i,j,k)=k2*vkr/kk-v2r
	    vi(i,j,k)=k2*vki/kk-v2i
	    wr(i,j,k)=k3*vkr/kk-v3r
	    wi(i,j,k)=k3*vki/kk-v3i
        else
         ur(i,j,k)=0.d0
         ui(i,j,k)=0.d0
         vr(i,j,k)=0.d0
         vi(i,j,k)=0.d0
         wr(i,j,k)=0.d0
         wi(i,j,k)=0.d0
        endif	  
   
       enddo
       enddo
       enddo
         
         call  FFT3D(nx,ny,nz,ur,ui,1)
         call  FFT3D(nx,ny,nz,vr,vi,1)
         call  FFT3D(nx,ny,nz,wr,wi,1)

         Er=0.d0
         Ei=0.d0
         do k=1,nz
         do j=1,ny
         do i=1,nx
          Er=Er+1.d0/2.d0*(ur(i,j,k)*ur(i,j,k)+vr(i,j,k)*vr(i,j,k)+  wr(i,j,k)*wr(i,j,k))
          Ei=Ei +1.d0/2.d0*(ui(i,j,k)*ui(i,j,k)+vi(i,j,k)*vi(i,j,k)+  wi(i,j,k)*wi(i,j,k))
         enddo
         enddo
         enddo

         Er=Er/(nx*ny*nz)
         Ei=Ei/(nx*ny*nz)
         print*, 'Er,Ei=',Er,Ei
         sr=sqrt(Ek0/Er)
         do k=1,nz
         do j=1,ny
         do i=1,nx
         ur(i,j,k)=ur(i,j,k)*sr
         vr(i,j,k)=vr(i,j,k)*sr
         wr(i,j,k)=wr(i,j,k)*sr
         enddo
         enddo
         enddo

!---------write init file (OpenCFD type) --------------------
!     ui is density and temperature
      do k=1,nz
      do j=1,ny
      do i=1,nx
        ui(i,j,k)=1.d0
       enddo
       enddo
       enddo
       
       print*, "write 3d file opencfd0.dat ..."
       open(99,file="opencfd0.dat",form="unformatted")
        write(99) 0, 0.d0
        call write3d(99,nx,ny,nz,ui)           ! density  (==1)
        call write3d(99,nx,ny,nz,ur)
        call write3d(99,nx,ny,nz,vr)
        call write3d(99,nx,ny,nz,wr)
        call write3d(99,nx,ny,nz,ui)          ! temperature (==1)
       close(99)
       
        deallocate(ur,ui,vr,vi,wr,wi)       
        end

!-----------------------------------------------------------------------
        subroutine write3d(fno,nx,ny,nz,f)
        implicit none
        integer:: nx,ny,nz,k,fno
        real*8:: f(nx,ny,nz)
        do k=1,nz
            write(fno) f(:,:,k)
        enddo
        end
        

!C    --------------- ---------------------------------------------------
!c               FFT FOR 1-D PROBLERM
!c      ****************************************************************
        SUBROUTINE FFT(N,FR,FI,SIGN)
        implicit double precision (a-h,o-z)
        dimension  FR(n),FI(n)
        integer sign
        MR=0
        NN=N-1
        DO 2 M=1,NN
         L=N
 1       L=L/2
        IF(MR+L.GT.NN) GO TO 1
        MR=MOD(MR,L)+L
        IF(MR.LE.M) GO TO 2
        TR=FR(M+1)
        FR(M+1)=FR(MR+1)
        FR(MR+1)=TR
        TI=FI(M+1)
        FI(M+1)=FI(MR+1)
        FI(MR+1)=TI
 2       CONTINUE
        L=1
 3       ISTEP=2*L
        EL=L
        DO 4 M=1,L
        A=3.1415926535897d0*dfloat(M-1)/EL
        WRo=COS(A)
        WIo=SIGN*SIN(A)
         DO 4 I=M,N,ISTEP
         J=I+L
         TR=WRo*FR(J)-WIo*FI(J)
         TI=WRo*FI(J)+WIo*FR(J)
         FR(J)=FR(I)-TR
         FI(J)=FI(I)-TI
         FR(I)=FR(I)+TR
 4        FI(I)=FI(I)+TI
        L=ISTEP
        IF(L.LT.N) GO TO 3
        IF(SIGN.GT.0.0) RETURN
        XN=1.d0/DFLOAT(N)
        DO 10 I=1,N
          FR(I)=FR(I)*XN
 10      FI(I)=FI(I)*XN
        RETURN
        END


!c----------------------------------------
       SUBROUTINE FFT2D(N1,N2,XR,XI,SIGN)
        implicit double precision (a-h,o-z)
        DIMENSION XR(n1,n2),XI(n1,n2)
        dimension FR(2048),FI(2048)
        integer sign

          DO  I=1,N1
         
          DO  J=1,N2
           FR(J)=XR(I,J)
           FI(J)=XI(I,J)
          ENDDO
          
         CALL FFT(N2,FR,FI,SIGN)
           DO  J=1,N2
           XR(I,J)=FR(J)
           XI(I,J)=FI(J)
          ENDDO
         ENDDO
         
         DO  J=1,N2
          DO  I=1,N1
          FR(I)=XR(I,J)
          FI(I)=XI(I,J)
          ENDDO
          CALL FFT(N1,FR,FI,SIGN)
          DO  I=1,N1
          XR(I,J)=FR(I)
          XI(I,J)=FI(I)
         ENDDO
         ENDDO
         
       RETURN
       END
!c =====================================
       SUBROUTINE FFT3D(N1,N2,N3,XR,XI,SIGN)
        implicit double precision (a-h,o-z)
        DIMENSION XR(n1,n2,n3),XI(n1,n2,n3)
        dimension FR(2048),FI(2048)
        integer sign
          do i=1,n1
          do j=1,n2
          do k=1,n3
           FR(k)=XR(i,j,k)
           FI(k)=XI(i,j,k)
          enddo
          call FFT(N3,FR,FI,SIGN)
          do  k=1,n3
          XR(i,j,k)=FR(k)
          XI(i,j,k)=FI(k)
          enddo
          enddo
          enddo
           do i=1,n1
          do k=1,n3
          do j=1,n2
           FR(j)=XR(i,j,k)
           FI(j)=XI(i,j,k)
          enddo
          call FFT(N2,FR,FI,SIGN)
          do  j=1,n2
          XR(i,j,k)=FR(j)
          XI(i,j,k)=FI(j)
          enddo
          enddo
          enddo

          do k=1,n3
          do j=1,n2
          do i=1,n1
           FR(i)=XR(i,j,k)
           FI(i)=XI(i,j,k)
          enddo
          call FFT(N1,FR,FI,SIGN)
          do  i=1,n1
          XR(i,j,k)=FR(i)
          XI(i,j,k)=FI(i)
          enddo
          enddo
          enddo

       RETURN
       END
  
