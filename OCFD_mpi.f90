!OpenCFD ver 2, CopyRight by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
!MPI Subroutines, such as computational domain partation, MPI message send and recv   

     subroutine partation3d_mpi
! Domain partation----------------------------------------------------------------------------
     Use flow_para
     implicit none
     integer k,ka,ierr
     integer npx1,npy1,npz1,npx2,npy2,npz2,my_mod1
!---------------------------------------------------------------------------------------------
	   call mpi_comm_size(MPI_COMM_WORLD,np_size,ierr)        
	   if(np_size .ne. npx0*npy0*npz0) then
         if(my_id.eq.0) print*, 'The Number of total Processes is not equal to npx0*npy0*npz0 !'
         stop
        endif
         
        npx=mod(my_id,npx0)
        npy=mod(my_id,npx0*npy0)/npx0
        npz=my_id/(npx0*npy0)
!------commonicaters-----------------------------------------------------------------------------
      
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz*npx0*npy0+npy*npx0, 0,MPI_COMM_X,ierr)   ! 1-D 
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz*npx0*npy0+npx, 0,MPI_COMM_Y,ierr)
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npy*npx0+npx, 0,MPI_COMM_Z,ierr)
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz, 0,MPI_COMM_XY,ierr)            ! 2-D 
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npy, 0,MPI_COMM_XZ,ierr)
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npx, 0,MPI_COMM_YZ,ierr)
      

!------------------------------------------------------------------------------------------------
!      均匀分配网格， 如果nx_global不能被npx0整除，将余下的网格点分到靠前的节点
!------------------------------------------------------------------------------------------------
      nx=nx_global/npx0
      ny=ny_global/npy0
      nz=nz_global/npz0
      if(npx .lt. mod(nx_global,npx0)) nx=nx+1
      if(npy .lt. mod(ny_global,npy0)) ny=ny+1
      if(npz .lt. mod(nz_global,npz0)) nz=nz+1

!------npx=k的节点上x方向网格点的个数，起始位置 -----------------------------
       allocate(i_offset(0:npx0-1),i_nn(0:npx0-1),  &
	            j_offset(0:npy0-1),j_nn(0:npy0-1),  &
				k_offset(0:npz0-1),k_nn(0:npz0-1))
				
!--------------------------------------------------------------------   
     do k=0,npx0-1
        ka=min(k,mod(nx_global,npx0))
        i_offset(k)=int(nx_global/npx0)*k+ka+1
        i_nn(k)=nx_global/npx0
        if(k .lt. mod(nx_global,npx0)) i_nn(k)=i_nn(k)+1
     enddo

     do k=0,npy0-1
        ka=min(k,mod(ny_global,npy0))
        j_offset(k)=int(ny_global/npy0)*k+ka+1
        j_nn(k)=ny_global/npy0
        if(k .lt. mod(ny_global,npy0)) j_nn(k)=j_nn(k)+1
     enddo

     do k=0,npz0-1
        ka=min(k,mod(nz_global,npz0))
        k_offset(k)=int(nz_global/npz0)*k+ka+1
        k_nn(k)=nz_global/npz0
        if(k .lt. mod(nz_global,npz0)) k_nn(k)=k_nn(k)+1
     enddo
!--------------------------------------------------------------------------------
!-------New Data TYPE------------------------------------------------------------
       call New_MPI_datatype

!--------define proc id:  the right, lift, up, bottom, frint and backward  procs
       npx1=my_mod1(npx-1,npx0)
       npx2=my_mod1(npx+1,npx0)
       ID_XM1=npz*(npx0*npy0)+npy*npx0+npx1    ! -1 proc in x-direction
       ID_XP1=npz*(npx0*npy0)+npy*npx0+npx2    ! +1 proc in x-direction
       if(Para%Iperiodic_X .eq.0 .and. npx .eq. 0) ID_XM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npx0-1 node
       if(Para%Iperiodic_X .eq.0 .and. npx .eq. npx0-1) ID_XP1=MPI_PROC_NULL
      
       npy1=my_mod1(npy-1,npy0)
       npy2=my_mod1(npy+1,npy0)
       ID_YM1=npz*(npx0*npy0)+npy1*npx0+npx
       ID_YP1=npz*(npx0*npy0)+npy2*npx0+npx
       if(Para%Iperiodic_Y.eq.0 .and. npy .eq. 0) ID_YM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npy0-1 node
       if(Para%Iperiodic_Y.eq.0 .and. npy .eq. npy0-1) ID_YP1=MPI_PROC_NULL

       npz1=my_mod1(npz-1,npz0)
       npz2=my_mod1(npz+1,npz0)
       ID_ZM1=npz1*(npx0*npy0)+npy*npx0+npx
       ID_ZP1=npz2*(npx0*npy0)+npy*npx0+npx
       if(Para%Iperiodic_Z.eq.0 .and. npz .eq. 0) ID_ZM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npz0-1 node
       if(Para%Iperiodic_Z.eq.0 .and. npz .eq. npz0-1) ID_ZP1=MPI_PROC_NULL


!--------------------------------------------------------------


      end            

!--------------------------------------------------------------------------------
         function my_mod1(i,n)
         implicit none
         integer my_mod1,i,n
           if(i<0) then
             my_mod1=i+n
           else if (i>n-1) then
             my_mod1=i-n
           else
             my_mod1=i
           endif
         end
!-----------------------------------------------------------------------------------------------
! Send Recv non-continuous data using derivative data type
   subroutine New_MPI_datatype
   Use flow_para
   implicit none
   integer:: ierr,TYPE_tmp
  
   call MPI_TYPE_Vector(ny,LAP,nx+2*LAP,OCFD_DATA_TYPE,TYPE_LAPX1,ierr)
   call MPI_TYPE_Vector(LAP,nx,nx+2*LAP,OCFD_DATA_TYPE,TYPE_LAPY1,ierr)
   call MPI_TYPE_Vector(LAP,nx,(nx+2*LAP)*(ny+2*LAP),OCFD_DATA_TYPE,TYPE_LAPZ1,ierr)

   call MPI_TYPE_Vector(ny,nx,nx+2*LAP,OCFD_DATA_TYPE,TYPE_tmp,ierr)

   call MPI_TYPE_HVector(nz,1,(nx+2*LAP)*(ny+2*LAP)*OCFD_REAL_KIND,TYPE_LAPX1,TYPE_LAPX2,ierr)
   call MPI_TYPE_HVector(nz,1,(nx+2*LAP)*(ny+2*LAP)*OCFD_REAL_KIND,TYPE_LAPY1,TYPE_LAPY2,ierr)
   call MPI_TYPE_HVector(LAP,1,(nx+2*LAP)*(ny+2*LAP)*OCFD_REAL_KIND,TYPE_tmp,TYPE_LAPZ2,ierr)

   call MPI_TYPE_COMMIT(TYPE_LAPX1,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPY1,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPZ1,ierr)

   call MPI_TYPE_COMMIT(TYPE_LAPX2,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPY2,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPZ2,ierr)

   call MPI_barrier(MPI_COMM_WORLD,ierr)
   end
!-----------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!----Form a Global index, get the node information and local index 
   subroutine get_i_node(i_global,node_i,i_local)
     Use Para_mpi
     implicit none
     integer i_global,node_i,i_local,ia

     node_i=npx0-1
     do ia=0,npx0-2
      if(i_global.ge.i_offset(ia) .and. i_global .lt. i_offset(ia+1) ) node_i=ia
     enddo
      i_local=i_global-i_offset(node_i)+1
     end
!-------------------------------------------------------------------------------
    subroutine get_j_node(j_global,node_j,j_local)
      Use Para_mpi
      implicit none
      integer j_global,node_j,j_local,ja
       node_j=npy0-1
       do ja=0,npy0-2
        if(j_global.ge.j_offset(ja) .and. j_global .lt. j_offset(ja+1) ) node_j=ja
       enddo
         j_local=j_global-j_offset(node_j)+1
     end
!-----------------------------------------------------------------------------------
    subroutine get_k_node(k_global,node_k,k_local)
       Use Para_mpi
       implicit none
       integer k_global,node_k,k_local,ka

       node_k=npz0-1
       do ka=0,npz0-2
        if(k_global .ge. k_offset(ka) .and. k_global .lt. k_offset(ka+1) ) node_k=ka
       enddo
       k_local=k_global-k_offset(node_k)+1
    end

!------------------------------------------------------------------------------------
       function get_id(npx1,npy1,npz1)
       Use Para_mpi
       implicit none
       integer get_id,npx1,npy1,npz1
       get_id=npz1*(npx0*npy0)+npy1*npx0+npx1
       return
       end
!-------------------------------------------------------------------------------------
!  Message send and recv at inner boundary (or 'MPI boundary')
       subroutine exchange_boundary_xyz(f)
       Use flow_para
       implicit none
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
       call exchange_boundary_x(f)
       call exchange_boundary_y(f)
       call exchange_boundary_z(f)
       return
       end

!=========================================================================================================
!  Boundary message communication (exchange bounary message)  
!=========================================================================================================

! mpi message send and recv, using user defined data type
       subroutine exchange_boundary_x(f)
       Use flow_para
       implicit none
       integer Status(MPI_status_Size),ierr
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      call MPI_Sendrecv(f(1,1,1),1,TYPE_LAPX2,  ID_XM1, 9000,   &
          f(nx+1,1,1),1, TYPE_LAPX2, ID_XP1, 9000,MPI_COMM_WORLD,Status,ierr)   
      call MPI_Sendrecv(f(nx+1-LAP,1,1),1,TYPE_LAPX2,  ID_XP1, 8000,   &           
          f(1-LAP,1,1),1,TYPE_LAPX2,  ID_XM1, 8000,MPI_COMM_WORLD,Status,ierr)   

      end
!------------------------------------------------------
       subroutine exchange_boundary_y(f)
       Use flow_para
       implicit none
       integer Status(MPI_status_Size),ierr
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)   

     call MPI_Sendrecv(f(1,1,1),1, TYPE_LAPY2,  ID_YM1,9000,   &
          f(1,ny+1,1),1,TYPE_LAPY2, ID_YP1, 9000,MPI_COMM_WORLD,Status,ierr)
     call MPI_Sendrecv(f(1,ny+1-LAP,1),1, TYPE_LAPY2, ID_YP1, 8000,  &
          f(1,1-LAP,1),1, TYPE_LAPY2,  ID_YM1, 8000,MPI_COMM_WORLD,Status,ierr)
     end
!------------------------------------------------------------
       subroutine exchange_boundary_z(f)
       Use flow_para
       implicit none
       integer Status(MPI_status_Size),ierr
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
     
        call MPI_Sendrecv(f(1,1,1),1,TYPE_LAPZ2, ID_ZM1,9000,   &
             f(1,1,nz+1),1, TYPE_LAPZ2, ID_ZP1, 9000,MPI_COMM_WORLD,Status,ierr)
        call MPI_Sendrecv(f(1,1,nz+1-LAP),1,TYPE_LAPZ2, ID_ZP1, 8000,   &
             f(1,1,1-LAP),1,TYPE_LAPZ2, ID_ZM1, 8000,MPI_COMM_WORLD,Status,ierr)
        end

!-----------------------------------------------------------------------

