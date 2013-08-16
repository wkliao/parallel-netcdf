!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$

! This program tests large files (> 4 GB) in PnetCDF. 

program tst_io
  use mpi
  use pnetcdf
  implicit none
  integer(kind=8), parameter :: prsz1 = 50, prsz2 = 50, &
       prsz3 = 50, prsz4 = 50, repct = 10
  integer :: i1, i2, i3, i4, j, k, ticksPerSec
  real :: psr
  integer :: clockRate, err
  integer :: start, now, wrint1, ncint1, wrint2, ncint2, &
       wrint3, ncint3, iosb, iosn, size
  real, dimension (prsz1, prsz2, prsz3, prsz4) :: s, t, u, v, w, x, y, z
  character(len = *), parameter :: nclFilenm1 = 'tst_io1.nc', &
       nclFilenm2 = 'tst_io2.nc', nclFilenm3 = 'tst_io3.nc', &
       nclFilenm4 = 'tst_io4.nc', nclFilenm5 = 'tst_io5.nc', &
       nclFilenm6 = 'tst_io6.nc', nclFilenm7 = 'tst_io7.nc', &
       nclFilenm8 = 'tst_io8.nc', nclFilenm9 = 'tst_io9.nc', &
       nclFilenm10 = 'tst_io10.nc', nclFilenm11 = 'tst_io11.nc'
  ! needed for netcdf
  integer :: ncid, x1id, x2id, x3id, x4id, vrid
  integer :: vrids, vridt, vridu, vridv, vridw, vridx, vridy, vridz

  psr = 1.7/real(prsz1)
  write(*,"(A)",advance="no") '*** TESTING Fortran tst_io.f90 '
  call MPI_Init(err)
  ! print *, "Starting data initialization."
  size = (prsz1 * prsz2 * prsz3 * prsz4 )/ 250000
  do i1 = 1, prsz1
     do i2 = 1, prsz2
        do i3 = 1, prsz3 ! Jackson Pollock it is not
           do i4 = 1, prsz4
              x(i1, i2, i3, i4) = sin(i1*psr)*(0.5 + cos(i2*psr))+(psr/i3)+ i4/(10.0*prsz4)
              s(i1, i2, i3, i4) = x(i1, i2, i3, i4) - 5.0
              t(i1, i2, i3, i4) = x(i1, i2, i3, i4) - 4.0
              u(i1, i2, i3, i4) = x(i1, i2, i3, i4) - 3.0
              v(i1, i2, i3, i4) = x(i1, i2, i3, i4) - 2.0
              w(i1, i2, i3, i4) = x(i1, i2, i3, i4) - 1.0
              y(i1, i2, i3, i4) = x(i1, i2, i3, i4) + 1.0
              z(i1, i2, i3, i4) = x(i1, i2, i3, i4) + 2.0
           enddo
        enddo
     enddo
  enddo

  call setupNetCDF (nclFilenm1, ncid, vrid, x, prsz1, prsz2, prsz3, prsz4, &
       x1id, x2id, x3id, x4id, NF90_CLOBBER, 20)
  call system_clock(start)
  call check(nfmpi_begin_indep_data(ncid), 11)
  call check (NF90MPI_PUT_VAR(ncid, vrid, x), 18)
  call check(nfmpi_end_indep_data(ncid), 11)
  call system_clock(now)
  ncint1 = now - start
!   print 3, size, "MB"," netcdf write = ", ncint1 * clockRate, &
!        real(ncint1)/real (wrint1)
3 format("Time for", i5, a, a25, i7, " msec. Spd ratio = ", f5.2)

  call check (NF90MPI_CLOSE(ncid), 14)

  call system_clock(start)
  do i1 = 1, repct
     call setupNetCDF (nclFilenm1, ncid, vrid, x, prsz1, prsz2, prsz3, prsz4, &
          x1id, x2id, x3id, x4id, NF90_CLOBBER, 130)
     call check(nfmpi_begin_indep_data(ncid), 11)
     call check (NF90MPI_PUT_VAR(ncid, vrid, x), 23 + i1)
     call check(nfmpi_end_indep_data(ncid), 11)
     call check (NF90MPI_CLOSE(ncid), 15)
  enddo
  call system_clock(now)
  ncint2 = now - start
!   print 4, size, repct, " repeated netcdf writes = ", ncint2 * clockRate, &
!        real(ncint2)/real(wrint2);
4 format("Time for", i5, "MB", i3, a22, i7, " msec. Spd ratio = ", f5.2)

!  call system_clock(start)
!  call setupNetCDF (nclFilenm3, ncid, vrids, s, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 20)
!  call setupNetCDF (nclFilenm4, ncid, vridt, t, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 30)
!  call setupNetCDF (nclFilenm5, ncid, vridu, u, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 40)
!  call setupNetCDF (nclFilenm6, ncid, vridv, v, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 50)
!  call setupNetCDF (nclFilenm7, ncid, vridw, w, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 60)
!  call setupNetCDF (nclFilenm8, ncid, vridx, x, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 70)
!  call setupNetCDF (nclFilenm9, ncid, vridy, y, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 80)
!  call setupNetCDF (nclFilenm10, ncid, vridz, z, prsz1, prsz2, prsz3, prsz4, &
!       x1id, x2id, x3id, x4id, NF90_CLOBBER, 90)
!  call check(nfmpi_begin_indep_data(ncid), 11)
!  call check (NF90MPI_PUT_VAR(ncid, vrids, s), 118)
!  call check (NF90MPI_PUT_VAR(ncid, vridt, t), 119)
!  call check (NF90MPI_PUT_VAR(ncid, vridu, u), 120)
!  call check (NF90MPI_PUT_VAR(ncid, vridv, v), 121)
!  call check (NF90MPI_PUT_VAR(ncid, vridw, w), 122)
!  call check (NF90MPI_PUT_VAR(ncid, vridx, x), 123)
!  call check (NF90MPI_PUT_VAR(ncid, vridy, y), 124)
!  call check (NF90MPI_PUT_VAR(ncid, vridz, z), 125)
!  call check(nfmpi_end_indep_data(ncid), 11)
!  call system_clock(now)
!  ncint3 = now - start
!  call check (NF90MPI_CLOSE(ncid), 16)
!   print 4, size, 8, " netcdf file writes = ", ncint3 * clockRate, &
!        real(ncint3)/real(wrint3);

  write(*,"(A)") '                                    ------ pass'

  call MPI_Finalize(err)

contains
  subroutine check (st, n) ! checks the return error code
    integer, intent (in) :: st, n
    if ((n < 10.and.st /= 0).or.(n > 10.and.st /= NF90_noerr))then
       print *, "I/O error at", n, " status = ", st
       stop 2
    endif
  end subroutine check

  subroutine setupNetCDF(fn, nc, vr, vrnam, d1, d2, d3, d4, do1, do2, &
       do3, do4, stat, deb)
    character(len = *), intent(in) :: fn
    real, dimension (:, :, :, :), intent (in) :: vrnam

    integer(kind=8), intent(in) :: d1, d2, d3, d4
    integer, intent(in) :: stat, deb
    integer, intent(out) :: do1, do2, do3, do4, vr
    integer, intent(inout) :: nc
    integer, dimension(4) :: dimids (4)

    call check (NF90MPI_CREATE (MPI_COMM_WORLD, fn, stat, MPI_INFO_NULL, nc), deb + 1)
    call check (NF90MPI_DEF_DIM(nc, "d1", d1, do1), deb + 2)
    call check (NF90MPI_DEF_DIM(nc, "d2", d2, do2), deb + 3)
    call check (NF90MPI_DEF_DIM(nc, "d3", d3, do3), deb + 4)
    call check (NF90MPI_DEF_DIM(nc, "d4", d4, do4), deb + 5)

    dimids = (/ do1, do2, do3, do4 /)
    call check (NF90MPI_DEF_VAR(nc, "data", NF90_REAL, dimids, vr), deb + 6)
    call check (NF90MPI_ENDDEF (nc), deb + 7)

  end subroutine setupNetCDF

end program tst_io