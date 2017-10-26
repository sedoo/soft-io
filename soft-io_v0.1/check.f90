subroutine check(status)

  integer, intent ( in) :: status
  character*80 strdum

  if ((status /= nf_noerr).and.(status /= 0)) then
     !       strdum=NF_STRERROR(status)
     !       print *, trim(strdum)
     print*,status,nf_noerr
     stop "Stopped - Netcdf Error"
  end if

end subroutine check
