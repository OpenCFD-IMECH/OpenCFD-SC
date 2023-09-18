!---This Code can interpolate data in z-direction, which should be a periodic dimension-----------
!-----------------------------------------------------------------------------------------
	  implicit doubleprecision(a-h,o-z)
!----all variables with footmark "1" are variables in old mesh; with "2" are in new mesh
	  real*8,allocatable,dimension(:):: z1,z2 ,v1,v2,w1,w2   
      real*8,allocatable,dimension(:,:,:):: f1,f2 
      integer,allocatable,dimension(:) :: ks_k
      real*8,allocatable,dimension(:,:)::  Ak
	  integer,parameter:: Periodic=1,Non_Periodic=0
      print*, "Interpolation flow in Z- direction, opencfd.dat --> opencfd-new.dat"
	  print*, "please input nx,ny "
	  read(*,*) nx,ny
	  print*, "Please input nz1 (origent mesh), nz2 (new mesh)"
	  read(*,*) nz1, nz2
! ---------------------------------------------------------

	  allocate( z1(nz1),z2(nz2))  !Coordinates
	  allocate ( f1(nx,ny,nz1), f2(nx,ny,nz2) )
	  allocate ( ks_k(nz2) )
      allocate (Ak(6,nz2),w1(nz1),w2(nz2))
!------------------------------------------------------------   

	  do k=1,nz1
	   z1(k)=(k-1.d0)/nz1
	  enddo

	  do k=1,nz2
	   z2(k)=(k-1.d0)/nz2
	  enddo

!--------------------------------------------------------

      call get_ks_Ai(nz1,nz2,z1,z2,ks_k,Ak,Periodic,1.d0)   ! 1.d0 Z-periodic length
      
	  print*, ks_k
      print*         
	  print*, Ak


	  open(44,file="opencfd.dat",form="unformatted")
	  open(55,file="opencfd-new.dat",form="unformatted")
      read(44) Istep,tt
      print*, "Istep,tt=", Istep,tt
      write(55) Istep,tt

	  do m=1,5
        print*, "interpolation ... ",m
        call read3d(44,f1,nx,ny,nz1)

!---------------interpolation in z- direction
         do j=1,ny
         do i=1,nx
          w1=f1(i,j,1:nz1)           
	      call interpolation6(nz1,nz2,ks_k,Ak,w1,w2,Periodic)
	      f2(i,j,1:nz2)=w2
         enddo
		 enddo

!------------------------------------------          
 	     call write3d(55,f2,nx,ny,nz2)
 	   enddo
      close(55)


!------------test the flow -------------------------------------------
       i=nx/2
       open(44,file="flow2d-yz-1.dat")
       write(44,*) "variables=x,y,T"
       write(44,*) "zone i= ",ny ,  " j= ", nz1
	   do k=1,nz1
	   do j=1,ny
	   write(44,"(3f15.6)") (j-1.d0)/(ny-1.d0),z1(k),f1(i,j,k)
	   enddo
	   enddo
       close(44)
         open(44,file="flow2d-yz-2.dat")
         write(44,*) "variables=x,y,T"
         write(44,*) "zone i= ", ny,  " j= ", nz2
	     do k=1,nz2
	     do j=1,ny
	     write(44,"(3f15.6)") (j-1.d0)/(ny-1.d0),z2(k),f2(i,j,k)
	     enddo
	     enddo
        close(44)

!--------------------------------------------------
      deallocate(f1,f2,z1,z2,Ak,ks_k)
       end




!-------------------------------------------------------
     subroutine get_ks_Ai(nx1,nx2,xx1,xx2,ks_i,Ai,IF_Periodic,Z_periodic)
	 implicit doubleprecision(a-h,o-z)
	 integer nx1,nx2,ks_i(nx2)
	 real*8 xx1(nx1),xx2(nx2),Ai(6,nx2)
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
           Ai(k,i)=1.d0

!	       ik=MOD_N(ks_i(i)+k-3,nx1)
	 	   ik=ks_i(i)+k-3
           if(ik < 1) then
		     ik=ik+nx1
		     xxa=xx1(ik)-Z_periodic
           else if(ik > nx1) then
		     ik=ik-nx1
		     xxa=xx1(ik)+Z_periodic
		   else
             xxa=xx1(ik)
		   endif

	   do km=ka,kb
	      ikm=ks_i(i)+km-3
!	      ikm=MOD_N(ks_i(i)+km-3,nx1)
          if(ikm > nx1) then
		   ikm=ikm-nx1
		   xxb=xx1(ikm)+Z_periodic
          else if( ikm < 1) then
		   ikm=ikm+nx1
		   xxb=xx1(ikm)-Z_periodic
		  else
		   xxb=xx1(ikm)
          endif

	     if(km .ne. k) then
!	      Ai(k,i)=Ai(k,i)*(xx2(i)-xx1(ikm))/(xx1(ik)-xx1(ikm))
	      Ai(k,i)=Ai(k,i)*(xx2(i)-xxb)/(xxa-xxb)
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
!---------------------------------------
    function MOD_N(k,N)
	 integer k,N,MOD_N
	 MOD_N=k
	 if(k .lt. 1) MOD_N=k+N
	 if(K .gt. N) MOD_N=k-N
	end

