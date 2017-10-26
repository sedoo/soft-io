! $Id: ioerror.f,v 3.5 2000/02/08 15:16:39 bmy Rel $
SUBROUTINE IOERROR( ERROR_NUM, UNIT, ROUTINE_NAME )
  !
  !*****************************************************************************
  !  Subroutine IOERRROR (bmy, 5/28/99, 2/7/00) tests for I/O errors.  If an 
  !  error condition exists, IOERROR will print out a warning message and the 
  !  error number, and then will halt execution.
  !
  !  On SGI systems, the user can use the "explain" command to get a more
  !  in-depth explanation of the I/O error condition.
  !
  !  Arguments as input:
  !  ===========================================================================
  !  (1) ERROR_NUM    : I/O error number (output from the IOSTAT flag)
  !  (2) UNIT         : Unit # of the file where the I/O error occurred
  !  (3) ROUTINE_NAME : Name of the routine in which the error occurred
  !
  !  NOTES:
  !  (1) Now flush the standard output buffer before stopping.  
  !      Also updated comments. (bmy, 2/7/00)
  !*****************************************************************************
  !  
  IMPLICIT NONE

  !**** Arguments
  INTEGER,             INTENT(IN) :: ERROR_NUM, UNIT
  CHARACTER ( LEN=* ), INTENT(IN) :: ROUTINE_NAME 

  !**** Local variables
  CHARACTER ( LEN=70 )            :: LINE
  !
  !*****************************************************************************
  !  IOERROR begins here!
  !
  !  Print a warning message and flush the std output buffer, so that the 
  !  error message will show up if output is piped to a log file
  !*****************************************************************************
  !
  WRITE( 6, 90 ) 
90 FORMAT( '=========================================',&
          '=========================================' ) 

  WRITE( 6, 100 ) ERROR_NUM, UNIT, TRIM( ROUTINE_NAME )
100 FORMAT( 'I/O error number ', i8, ' in file unit ', i8, /, &
          'Encountered in routine ', a )

  WRITE( 6, 90 ) 

  CALL FLUSH( 6 )
  !
  !*****************************************************************************
  !  Halt execution
  !*****************************************************************************
  !
  STOP

END SUBROUTINE IOERROR
