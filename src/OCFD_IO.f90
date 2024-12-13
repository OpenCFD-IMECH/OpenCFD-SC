!Read & save file
!Copyright by Li Xinliang, lixl@lnm.imech.ac.cn
!-----------------------------------------------------------------------------------------
! read file, OpenCFD-1.5
  subroutine read_flow_data
     use flow_data
     implicit none
     character(len=100):: filename
     integer:: ierr
     logical :: file_exist
!-----------------------------------------------------------
     filename = "opencfd.dat"

     ! Check if the file exists
     inquire( file=trim(filename), exist=file_exist )

     if(.not. file_exist) then
      call flowfield_init('TGV')
      return
     endif
     
     if (my_id .eq. 0) then
        print *, 'read initial data file: ', filename
        open (55, file=filename, form='unformatted')
        read (55) Istep, tt
        print *, Istep, tt
     end if

     call MPI_bcast(Istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_bcast(tt, 1, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, ierr)
     call read_3d1(55, d)
     call read_3d1(55, u)
     call read_3d1(55, v)
     call read_3d1(55, w)
     call read_3d1(55, T)

     if (my_id .eq. 0) then
        close (55)
        print *, 'read data ok'
     end if
  end

!================================================================================
! save file, OpenCFD 1.5
  subroutine save_flow_data
     use flow_data
     implicit none
     character(len=100) filename1
     write (filename1, "('OCFD'I8.8'.dat')") Istep

     if (my_id .eq. 0) then
        print *, 'write data file:', filename1
        open (155, file=filename1, form='unformatted')
        write (155) Istep, tt
     end if
     call write_3d1(155, d)
     call write_3d1(155, u)
     call write_3d1(155, v)
     call write_3d1(155, w)
     call write_3d1(155, T)
     if (my_id .eq. 0) then
        close (155)
        print *, 'write data ok'
     end if
  end
!-------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
  subroutine write_3d1(file_no, U)
     use flow_para
     implicit none
     integer file_no, i, j, k
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND), allocatable:: U1(:, :, :)
     allocate (U1(nx, ny, nz))
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        U1(i, j, k) = U(i, j, k)
     end do
     end do
     end do
     call write_3d(file_no, U1, nx, ny, nz, nx_global, ny_global, nz_global)
     deallocate (U1)
  end

  subroutine read_3d1(file_no, U)
     use flow_para
     implicit none
     integer file_no, i, j, k
     real(kind=OCFD_REAL_KIND):: U(1 - LAP:nx + LAP, 1 - LAP:ny + LAP, 1 - LAP:nz + LAP)
     real(kind=OCFD_REAL_KIND), allocatable:: U1(:, :, :)
     allocate (U1(nx, ny, nz))
     call read_3d(file_no, U1, nx, ny, nz, nx_global, ny_global, nz_global)
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        U(i, j, k) = U1(i, j, k)
     end do
     end do
     end do
     deallocate (U1)
  end

