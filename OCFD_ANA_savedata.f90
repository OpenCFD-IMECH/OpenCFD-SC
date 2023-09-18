
!------ code for saving data -----------------------------------------------
  subroutine ana_savedata(ana_data)  
   Use flow_data     
   implicit none
   integer i,j,k,savetype
   real(kind=OCFD_REAL_KIND):: ana_data(50)
   integer,parameter:: Save_YZ=1,Save_XZ=2,Save_XY=3,Save_Point=4,Save_Block=5  

   savetype=nint(ana_data(1))
   select case(savetype)
    case(Save_YZ) 
	  call ana_saveplaneYZ(ana_data)
    case(Save_XZ)
	  call ana_saveplaneXZ(ana_data)
    case(Save_XY)
      call ana_saveplaneXY(ana_data)
    case(Save_Point)
      call ana_savePoints(ana_data)
    case(Save_Block)	  
	  call ana_saveblock(ana_data)
   end select
  end 
 
 
!------------------------------------------------------------------------------
! Save plane i_global=i0 to i0+bandwidth  
 subroutine ana_saveplaneYZ(ana_data)
  Use flow_data     
  implicit none
  integer i,j,k,ka,ia,ierr,points,bandwidth,i_point
  real(kind=OCFD_REAL_KIND):: ana_data(50)
  character(len=50) filename

!------------------------------------------------------------
     points=nint(ana_data(2))                 ! point number
     bandwidth=nint(ana_data(3))              ! bandwidth of the saving i- location 
     
    do ka=1,points

      if(my_id.eq.0) then
        print*, 'save data ....',Istep,tt,ka
        write(filename,'("savedata-YZ"I3.3".dat")') ka
        open(66,file=filename,form='unformatted',position='append')
        write(66) Istep,tt
      endif
     
	 i_point=nint(ana_data(3+ka))                  ! i-location of the saving points
     do ia=i_point,i_point+bandwidth-1
       call write_2d_YZa(66,ia,d)
       call write_2d_YZa(66,ia,u)
       call write_2d_YZa(66,ia,v)
       call write_2d_YZa(66,ia,w)
       call write_2d_YZa(66,ia,T)
     enddo
     if(my_id.eq.0)    close(66)
    enddo
    end
!--------------------------------------------------------
 subroutine ana_saveplaneXZ(ana_data)
  Use flow_data     
  implicit none
  integer ja,ka,ierr,points,bandwidth,j_point
  real(kind=OCFD_REAL_KIND):: ana_data(50)
  character(len=50) filename

!------------------------------------------------------------
     points=nint(ana_data(2))
     bandwidth=nint(ana_data(3))
   
    do ka=1,points
      if(my_id.eq.0) then
        print*, 'save data ....',Istep,tt,ka
        write(filename,'("savedata-XZ"I3.3".dat")') ka
        open(66,file=filename,form='unformatted',position='append')
        write(66) Istep,tt
      endif
	  
     j_point=nint(ana_data(3+ka))
     do ja=j_point,j_point+bandwidth-1
       call write_2d_XZa(66,ja,d)
       call write_2d_XZa(66,ja,u)
       call write_2d_XZa(66,ja,v)
       call write_2d_XZa(66,ja,w)
       call write_2d_XZa(66,ja,T)
     enddo
    
       if(my_id.eq.0)    close(66)
      
     enddo

   end

!-----------------------------------------------------------
 subroutine ana_saveplaneXY(ana_data)
  Use flow_data     
  implicit none
  integer i,j,k,ka,ia,ierr,points,bandwidth,k_point
  real(kind=OCFD_REAL_KIND):: ana_data(50)
  character(len=50) filename

!------------------------------------------------------------
    points=nint(ana_data(2))
    bandwidth=nint(ana_data(3))
	 
    do ka=1,points

      if(my_id.eq.0) then
        print*, 'save data ....',Istep,tt,ka
        write(filename,'("savedata-XY"I3.3".dat")') ka
        open(66,file=filename,form='unformatted',position='append')
        write(66) Istep,tt
      endif
     k_point=nint(ana_data(3+ka))
     do ia=k_point,k_point+bandwidth-1
       call write_2d_XYa(66,ia,d)
       call write_2d_XYa(66,ia,u)
       call write_2d_XYa(66,ia,v)
       call write_2d_XYa(66,ia,w)
       call write_2d_XYa(66,ia,T)
     enddo
     if(my_id.eq.0)   close(66)
    enddo
    end

!--------------------------------------------------------
 subroutine ana_savePoints(ana_data)
  Use flow_data     
  implicit none
  integer k,ierr
  integer,save:: Kflag1=0,mpoints
  integer,save,dimension(10000)::ia,ja,ka
  real(kind=OCFD_REAL_KIND):: ana_data(50)

!------------------------------------------------------------
   if(Kflag1.eq.0) then
      Kflag1=1
         if(my_id.eq.0) then
          open(99,file="save_points_locations.dat")
          read(99,*) mpoints
          do k=1,mpoints
            read(99,*) ia(k),ja(k),ka(k)
          enddo
          close(99)
         endif
      call MPI_bcast(mpoints,1,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)
      call MPI_bcast(ia,10000,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)
      call MPI_bcast(ja,10000,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)
      call MPI_bcast(ka,10000,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)
   endif
 
      if(my_id.eq.0) then
        print*, 'save points ......',Istep,tt,mpoints
        open(66,file="save_points.dat",form='unformatted',position='append')
        write(66) tt
      endif
     
       call write_points(66,d,mpoints,ia,ja,ka)
       call write_points(66,u,mpoints,ia,ja,ka)
       call write_points(66,v,mpoints,ia,ja,ka)
       call write_points(66,w,mpoints,ia,ja,ka)
       call write_points(66,T,mpoints,ia,ja,ka)

       if(my_id.eq.0) then
       close(66)
       endif
   end
!--------------------------------------------------------
! save block data 
!-----------------------------------------------------------
 subroutine ana_saveblock(ana_data)
  Use flow_data     
  implicit none
  integer ib,ie,jb,je,kb,ke,filenumber_K
  character(len=50) filename
  real(kind=OCFD_REAL_KIND):: ana_data(50)
	
	 ib=nint(ana_data(2)) ;  ie=nint(ana_data(3))
	 jb=nint(ana_data(4)) ;  je=nint(ana_data(5))
	 kb=nint(ana_data(6)) ;  ke=nint(ana_data(7))
     filenumber_K=nint(ana_data(8))

     if(my_id.eq.0) then
     write(filename,'("savedata-block"I3.3".dat")') filenumber_K
     print*, 'save block data ....',filename, Istep,tt
     open(66,file=filename,form='unformatted',position='append')
     write(66) Istep,tt
     endif
     
	 call write_blockdata(66,d,ib,ie,jb,je,kb,ke)
     call write_blockdata(66,u,ib,ie,jb,je,kb,ke)
     call write_blockdata(66,v,ib,ie,jb,je,kb,ke)
     call write_blockdata(66,w,ib,ie,jb,je,kb,ke)
     call write_blockdata(66,T,ib,ie,jb,je,kb,ke)
	 if(my_id .eq. 0)  close(66)
end
  

    
!----------------