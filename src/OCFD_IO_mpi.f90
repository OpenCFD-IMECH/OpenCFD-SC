! Read and write files by MPI
! Copyright by Li Xinliang
!------------------------------------------------------------

!--------------------------------------------------------------------------
!  write f(1:nx_global, 1:ny_global, ka)   (f=d,u,v,w,T)
  subroutine write_2d_XYa(file_no, ka, U)
! write a xy-plane data, i.e U_globle(1:nx_global,1:ny_lobal,ka)
! Use MPI_REDUCE
     use flow_para
     implicit none
     integer file_no, ka, node_k, k_local, i, j, i1, j1, ierr
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :) ::  U2d, U0

     allocate (U2d(nx_global, ny_global), U0(nx_global, ny_global))
     U2d = 0.d0
     U0 = 0.d0
     call get_k_node(ka, node_k, k_local)
     if (npz .eq. node_k) then
        do j = 1, ny
           do i = 1, nx
              i1 = i_offset(npx) + i - 1
              j1 = j_offset(npy) + j - 1
              U2d(i1, j1) = U(i, j, k_local)
           end do
        end do
     end if
     call MPI_Reduce(U2d, U0, nx_global*ny_global, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     if (my_id .eq. 0) then
        write (file_no) U0
     end if
     deallocate (U2d, U0)
  end

!--------------------------------------------------------------
!-----Write a 2D Y-Z (j-k) plane from 3-D array
  subroutine write_2d_YZa(file_no, ia, U)
     use flow_para
     implicit none
     integer ierr, ia, node_i, i_local
     integer file_no, j, k, j1, k1
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :):: U2d, U0

!--------------------------------------------
     allocate (U2d(ny_global, nz_global), U0(ny_global, nz_global))
     U2d = 0.d0
     U0 = 0.d0

     call get_i_node(ia, node_i, i_local)
     if (npx .eq. node_i) then
        do k = 1, nz
        do j = 1, ny
           j1 = j_offset(npy) + j - 1
           k1 = k_offset(npz) + k - 1
           U2d(j1, k1) = U(i_local, j, k)
        end do
        end do
     end if

     call MPI_Reduce(U2d, U0, ny_global*nz_global, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     if (my_id .eq. 0) then
        write (file_no) U0
     end if
     deallocate (U2d, U0)
  end

!-------------------------------------------------
!----Write a 2d xz-plane from 3d array------------------------

  subroutine write_2d_XZa(file_no, ja, U)
     use flow_para
     implicit none
     integer file_no, ierr, i, k, i1, k1, ja, node_j, j_local
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :):: U2d, U0
     allocate (U2d(nx_global, nz_global), U0(nx_global, nz_global))

     U2d = 0.d0
     U0 = 0.d0
     call get_j_node(ja, node_j, j_local)
     if (npy .eq. node_j) then
        do k = 1, nz
           do i = 1, nx
              i1 = i_offset(npx) + i - 1
              k1 = k_offset(npz) + k - 1
              U2d(i1, k1) = U(i, j_local, k)
           end do
        end do
     end if
     call MPI_Reduce(U2d, U0, nx_global*nz_global, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     if (my_id .eq. 0) then
        write (file_no) U0
     end if
     deallocate (U2d, U0)
  end

!--------------------------------------------------

!----Write points from 3d array------------------------

  subroutine write_points(file_no, U, mpoints, ia, ja, ka)
     use flow_para
     implicit none
     integer mpoints, m, nrecv
     integer Status(MPI_status_Size), ierr
     integer file_no, node_i, node_j, node_k, i_local, j_local, k_local
     integer, dimension(mpoints):: ia, ja, ka
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND):: U1(mpoints)
!--------------------------------
     do m = 1, mpoints
        call get_i_node(ia(m), node_i, i_local)
        call get_j_node(ja(m), node_j, j_local)
        call get_k_node(ka(m), node_k, k_local)
        if (npx .eq. node_i .and. npy .eq. node_j .and. npz .eq. node_k) then
           call MPI_Bsend(U(i_local, j_local, k_local), 1, OCFD_DATA_TYPE, 0, 0, MPI_COMM_WORLD, ierr)
        end if

        if (my_id .eq. 0) then
           nrecv = node_k*npx0*npy0 + node_j*npx0 + node_i
           call MPI_Recv(U1(m), 1, OCFD_DATA_TYPE, nrecv, 0, MPI_COMM_WORLD, status, ierr)
        end if
     end do
     if (my_id .eq. 0) then
        write (file_no) U1
     end if
  end

!--------------------------------------------------
  subroutine read_3d(file_no, U, nx, ny, nz, nx0, ny0, nz0)
     use OCFD_precision
     use Para_mpi
     implicit none
     integer nx, ny, nz, nx0, ny0, nz0
     integer Status(MPI_status_Size), ierr, i1, j1, k1, nk, i0, j0, ia, i, j, k, kk
     integer file_no, npx1, npy1, npz1, recvcount
     real(kind=OCFD_REAL_KIND):: U(nx, ny, nz)
     real(kind=OCFD_REAL_KIND), allocatable:: buff2d(:, :), buff1(:, :), buff2(:), buff_recv(:)
     integer, allocatable:: sendcounts1(:), displs1(:), sendcounts2(:), displs2(:)
!---------------------------------------------------------------
     allocate (buff2d(nx0, ny0), buff1(nx0, ny), buff2(nx0*ny), buff_recv(nx*ny))
     allocate (sendcounts1(npy0), displs1(npy0), sendcounts2(npx0), displs2(npx0))

     if (my_id .eq. 0) print *, 'read 3d data ...'

     do npy1 = 0, npy0 - 1
        sendcounts1(npy1 + 1) = nx0*j_nn(npy1)
        if (npy1 .eq. 0) then
           displs1(npy1 + 1) = 0
        else
           displs1(npy1 + 1) = displs1(npy1) + sendcounts1(npy1)
        end if
     end do

     do npx1 = 0, npx0 - 1
        sendcounts2(npx1 + 1) = i_nn(npx1)*ny
        if (npx1 .eq. 0) then
           displs2(npx1 + 1) = 0
        else
           displs2(npx1 + 1) = displs2(npx1) + sendcounts2(npx1)
        end if
     end do

     do kk = 1, nz0
        call get_k_node(kk, nk, k1)

        if (my_id .eq. 0) then
           read (file_no) buff2d
        end if

        if (nk .ne. 0) then
           if (my_id .eq. 0) call MPI_Send(buff2d, nx0*ny0, OCFD_DATA_TYPE, nk*(npx0*npy0), 6666, MPI_COMM_WORLD, ierr)
           if (npx .eq. 0 .and. npy .eq. 0 .and. npz .eq. nk) &
              call MPI_recv(buff2d, nx0*ny0, OCFD_DATA_TYPE, 0, 6666, MPI_COMM_WORLD, status, ierr)
        end if

        if (npz .eq. nk) then
           if (npx .eq. 0) then
              call MPI_scatterv(buff2d, sendcounts1, displs1, OCFD_DATA_TYPE, buff1, nx0*ny, OCFD_DATA_TYPE, 0, MPI_COMM_Y, ierr)
           end if
           ! re-arrage data ......, transfer buff1 to buff2
           ia = 0
           do npx1 = 0, npx0 - 1
           do i = 1, i_nn(npx1)*ny
              i0 = i_offset(npx1) + mod(i - 1, i_nn(npx1))
              j0 = int((i - 1)/i_nn(npx1)) + 1
              buff2(ia + i) = buff1(i0, j0)
           end do
           ia = ia + i_nn(npx1)*ny
           end do

           call MPI_scatterv(buff2, sendcounts2, displs2, OCFD_DATA_TYPE, buff_recv, nx*ny, OCFD_DATA_TYPE, 0, MPI_COMM_X, ierr)

           do j1 = 1, ny
           do i1 = 1, nx
              U(i1, j1, k1) = buff_recv(i1 + (j1 - 1)*nx)
           end do
           end do
        end if
     end do
     deallocate (buff2d, buff1, buff2, buff_recv)
     deallocate (sendcounts1, displs1, sendcounts2, displs2)
  end

!------------------------------------------------------------------------------------------
  subroutine write_3d(file_no, U, nx, ny, nz, nx0, ny0, nz0)
     use OCFD_precision
     use Para_mpi
     implicit none
     integer nx, ny, nz, nx0, ny0, nz0
     integer Status(MPI_status_Size), ierr, i1, j1, k1, nk, i0, j0, ia, i, j, k, kk
     integer file_no, npx1, npy1, npz1, recvcount
     real(kind=OCFD_REAL_KIND):: U(nx, ny, nz)
     real(kind=OCFD_REAL_KIND), allocatable:: buff2d(:, :), buff1(:, :), buff2(:), buff_send(:)
     integer, allocatable:: recvcounts1(:), displs1(:), recvcounts2(:), displs2(:)
!---------------------------------------------------------------
     allocate (buff2d(nx0, ny0), buff1(nx0, ny), buff2(nx0*ny), buff_send(nx*ny))
     allocate (recvcounts1(npy0), displs1(npy0), recvcounts2(npx0), displs2(npx0))

     if (my_id .eq. 0) print *, 'write 3d data ...'

     do npy1 = 0, npy0 - 1
        recvcounts1(npy1 + 1) = nx0*j_nn(npy1)
        if (npy1 .eq. 0) then
           displs1(npy1 + 1) = 0
        else
           displs1(npy1 + 1) = displs1(npy1) + recvcounts1(npy1)
        end if
     end do

     do npx1 = 0, npx0 - 1
        recvcounts2(npx1 + 1) = i_nn(npx1)*ny
        if (npx1 .eq. 0) then
           displs2(npx1 + 1) = 0
        else
           displs2(npx1 + 1) = displs2(npx1) + recvcounts2(npx1)
        end if
     end do

     do kk = 1, nz0
        call get_k_node(kk, nk, k1)

        if (npz .eq. nk) then
           do j1 = 1, ny
           do i1 = 1, nx
              buff_send(i1 + (j1 - 1)*nx) = U(i1, j1, k1)
           end do
           end do
           call MPI_gatherv(buff_send, nx*ny, OCFD_DATA_TYPE, buff2, recvcounts2, displs2, OCFD_DATA_TYPE, 0, MPI_COMM_X, ierr)

           ia = 0
           do npx1 = 0, npx0 - 1
           do i = 1, i_nn(npx1)*ny
              i0 = i_offset(npx1) + mod(i - 1, i_nn(npx1))     ! i_offset(npx1)+(mod(i-1,i_nn(npx1))+1) -1
              j0 = int((i - 1)/i_nn(npx1)) + 1
              buff1(i0, j0) = buff2(ia + i)
           end do
           ia = ia + i_nn(npx1)*ny
           end do

           if (npx .eq. 0) then
              call MPI_gatherv(buff1, nx0*ny, OCFD_DATA_TYPE, buff2d, recvcounts1, displs1, OCFD_DATA_TYPE, 0, MPI_COMM_Y, ierr)
           end if
        end if

        if (nk .ne. 0) then
           if (npx .eq. 0 .and. npy .eq. 0 .and. npz .eq. nk) &
              call MPI_send(buff2d, nx0*ny0, OCFD_DATA_TYPE, 0, 6666, MPI_COMM_WORLD, ierr)
           if (my_id .eq. 0) call MPI_recv(buff2d, nx0*ny0, OCFD_DATA_TYPE, nk*(npx0*npy0), 6666, MPI_COMM_WORLD, status, ierr)
        end if

        if (my_id .eq. 0) then
           write (file_no) buff2d
        end if
     end do

     deallocate (buff2d, buff1, buff2, buff_send)
     deallocate (recvcounts1, displs1, recvcounts2, displs2)

  end

!--------------------------------Write blockdata from 3d array-------------------------------------------------

  subroutine write_blockdata(file_no, U, ib, ie, jb, je, kb, ke)
     use flow_para
     implicit none
     integer Status(MPI_status_Size), ierr
     integer::  file_no, ib, ie, jb, je, kb, ke, nx1, ny1, nz1, i, j, k, i0, j0, k0, i1, j1, k1
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND), allocatable, dimension(:, :, :):: U1, U0
!--------------------------------
     nx1 = ie - ib + 1; ny1 = je - jb + 1; nz1 = ke - kb + 1
     allocate (U1(nx1, ny1, nz1), U0(nx1, ny1, nz1))
     do k = 1, nz1
     do j = 1, ny1
     do i = 1, nx1
        U1(i, j, k) = 0.d0
        U0(i, j, k) = 0.d0
     end do
     end do
     end do

     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        k0 = k_offset(npz) + k - 1
        j0 = j_offset(npy) + j - 1
        i0 = i_offset(npx) + i - 1
        if (i0 >= ib .and. i0 <= ie .and. j0 >= jb .and. j0 <= je .and. k0 >= kb .and. k0 <= ke) then
           i1 = i0 - ib + 1
           j1 = j0 - jb + 1
           k1 = k0 - kb + 1
           U1(i1, j1, k1) = U(i, j, k)
        end if
     end do
     end do
     end do

     call MPI_Reduce(U1, U0, nx1*ny1*nz1, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

     if (my_id .eq. 0) then
        write (file_no) U0
     end if
     deallocate (U1, U0)
  end

!--------------------------------------------------

