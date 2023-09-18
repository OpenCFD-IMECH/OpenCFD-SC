! call FFTW for FFT
! ver 1.1 2011-9-25
! get Kinetic Energy Spectrum  using FFT3d 
     implicit none
!     include 'fftw_f77.i'  ! FFTW2.5
!     include "fftw3.f"     ! since FFTW3.0 the include file became "fftw3.f"
     integer,parameter::nx=128,ny=nx,nz=nx,Kmax=2*nx
     real*8,parameter:: PI=3.1415926535d0
     doublecomplex CI,out
     real*8,allocatable:: uu(:,:,:),ur(:,:,:),ui(:,:,:),Er(:,:,:)
     real*8 Ek(0:Kmax),Ea,Ek1,Ek2,Ek3,Ek4,tt,urms,uxrms,le,tmp,Re
     integer  plan,ierr,nn(3),i,j,k,k1,k2,k3,kk,Istep,i1
     CI=(0.d0,1.d0)

          nn(1)=nx
	  nn(2)=ny
	  nn(3)=nz

      allocate(ur(nx,ny,nz),ui(nx,ny,nz),Er(nx,ny,nz))
      allocate(uu(nx,ny,nz))
!--------------------------------------------------
      print*, "Comput 1-d Energy spectrum E(k) (shell-integrated energy spectrum)"
      print*, "Comput integrate length le, eddy turn-over time tao"
      print*, "please input Re ..."
      read(*,*) Re

      urms=0.d0
!      call fftwnd_f77_create_plan(plan,3,nn,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
      open(33,file='opencfd.dat',form='unformatted')
      read(33) Istep,tt
      print*, 'Istep=',Istep,'tt=',tt
      call read3d(33,nx,ny,nz,uu)  ! d
      call read3d(33,nx,ny,nz,uu)  ! u

       do k=1,nz
       do j=1,ny
       do i=1,nx
        urms=urms+uu(i,j,k)**2
       enddo
       enddo
       enddo

 !     call fftwnd_f77_one(plan,uf,out)   ! spectrum of u
      ur=uu 
      ui=0.d0
      call FFT3D(nx,ny,nz,ur,ui,-1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
       Er(i,j,k)=(ur(i,j,k)**2+ui(i,j,k)**2)/2.d0
      enddo
      enddo
      enddo

!-----get RMS of ux -------------------
    do  k=1,Nz
!         k1=k-1 ;    if(k1 .ge. Nz/2 ) k1=k1-Nz
    do  j=1,Ny
!         j1=j-1 ;    if(j1 .ge. Ny/2 ) j1=j1-Ny
    do  i=1,Nx
         i1=i-1 ;    if(i1 .ge. Nx/2 ) i1=i1-Nx
       tmp=ur(i,j,k)
       ur(i,j,k)=-ui(i,j,k)*i1
       ui(i,j,k)=tmp*i1
    enddo
    enddo
    enddo
   
    call FFT3D(nx,ny,nz,ur,ui,1)
     uxrms=0.d0
       do k=1,nz
       do j=1,ny
       do i=1,nx
        uxrms=uxrms+ur(i,j,k)**2
       enddo
       enddo
       enddo
       uxrms=sqrt(uxrms/(nx*ny*nz))
!-------------------------------------
      call read3d(33,nx,ny,nz,uu)  ! v
       do k=1,nz
       do j=1,ny
       do i=1,nx
        urms=urms+uu(i,j,k)**2
       enddo
       enddo
       enddo
      ur=uu
      ui=0.d0


      call FFT3D(nx,ny,nz,ur,ui,-1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
       Er(i,j,k)=Er(i,j,k)+(ur(i,j,k)**2+ui(i,j,k)**2)/2.d0
      enddo
      enddo
      enddo

 !     call fftwnd_f77_one(plan,vf,out)   ! spectrum of u
      call read3d(33,nx,ny,nz,uu)  ! w

       do k=1,nz
       do j=1,ny
       do i=1,nx
        urms=urms+uu(i,j,k)**2
       enddo
       enddo
       enddo
  
      ur=uu
      ui=0.d0
      call FFT3D(nx,ny,nz,ur,ui,-1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
       Er(i,j,k)=Er(i,j,k)+(ur(i,j,k)**2+ui(i,j,k)**2)/2.d0
      enddo
      enddo
      enddo

 !     call fftwnd_f77_one(plan,wf,out)   ! spectrum of u

! Ek(k) is averaged Kinetic Energy Spectrum in the interval [k-0.5,k+0.5) (k=sqrt(k1**2+k2**2+k3**2))
!      uf=uf/(1.d0*nx*ny*nz)  ! FFTW do not times 1.d0/(nx*ny*nz)
!      vf=vf/(1.d0*nx*ny*nz)
!      wf=wf/(1.d0*nx*ny*nz)
!--------------------------------
      urms=sqrt(urms/(3.d0*nx*ny*nz))
      print*, "statistics rms =", urms
 
      Ek=0.d0
      do k=1,nz
	  do j=1,ny
	  do i=1,nx
        k3=k-1
        if(k3 .ge. nz/2) k3=k-1-nz
        k2=j-1
        if(k2 .ge. ny/2) k2=j-1-ny
        k1=i-1
        if(k1 .ge. nx/2) k1=i-1-nx

         kk=int(sqrt(dble(k1*k1+k2*k2+k3*k3))+0.5d0)
         Ek(kk)=Ek(kk)+Er(i,j,k)
       enddo
       enddo
       enddo

	   open(66,file='Energy-spectrum.dat')
           Ek1=0.d0
           Ek2=0.d0
	   do k=1,Kmax
           Ek1=Ek1+Ek(k)
           Ek2=Ek2+Ek(k)/k
 	   write(66,*) k, Ek(k)
	   enddo
	   close(66)
           urms=sqrt(Ek1*2.d0/3.d0)
           le= 3.d0/4.d0*PI*Ek2/Ek1
           print*, "urms=", urms, "uxrms=", uxrms
           print*, "Integrated length le=", le
           print*, "Taylor scale (lamda)=", urms/uxrms
           print*, "Eddy turnover time tao=", le/urms
           print*, "Re_lamda=", urms*(urms/uxrms)*Re
	   end

! -----------------------------------------------
       subroutine read3d(file_no,nx,ny,nz,u)
       implicit none
       integer file_no,nx,ny,nz,k
       real*8 u(nx,ny,nz),u2d(nx,ny)
       do k=1,nz
       read(file_no) u2d
       u(:,:,k)=u2d
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
          DO 4 I=1,N1
          DO 3 J=1,N2
           FR(J)=XR(I,J)
 3         FI(J)=XI(I,J)
         call FFT(N2,FR,FI,SIGN)
          DO 4 J=1,N2
          XR(I,J)=FR(J)
 4         XI(I,J)=FI(J)
         DO 6 J=1,N2
          DO 5 I=1,N1
          FR(I)=XR(I,J)
 5         FI(I)=XI(I,J)
         call FFT(N1,FR,FI,SIGN)
          DO 6 I=1,N1
          XR(I,J)=FR(I)
 6         XI(I,J)=FI(I)
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
  
