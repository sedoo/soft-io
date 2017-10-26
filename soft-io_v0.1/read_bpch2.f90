SUBROUTINE READ_BPCH2( FILENAME, CATEGORY_IN, TRACER_IN, TAU0_IN, &
     IX,JX,LX,ARRAY ) 
!
!*****************************************************************************
!  Subroutine READ_BPCH2 reads a binary punch file (v. 2.0) and extracts
!  a data block that matches the given category, tracer, and tau value.
!  (bmy, 12/10/99)
!
!  Arguments as Input:
!  ===========================================================================
!  (1  ) FILENAME    : (CHARACTER) String for input file name
!  (2  ) CATEGORY_IN : (CHARACTER) Category name for the desired data block
!  (3  ) TRACER_IN   : (INTEGER  ) Tracer number for which to extract data
!  (4  ) TAU0_IN     : (REAL*8   ) TAU value for which to extract data
!  (5-7) IX, JX, LX  : (INTEGER  ) Dimensions of ARRAY (see below) 
!
  !  Arguments as Output:
  !  ===========================================================================
  !  (8  ) ARRAY       : (REAL*4   ) Array to hold extracted data values
  !
  !  NOTES:
  !  (1) Trap all I/O errors with subroutine IOERROR.F.
  !*****************************************************************************
  !
  IMPLICIT NONE

  include "CMN_SIZE"

  !**** Arguments
  INTEGER,           INTENT(IN)  :: IX, JX, LX, TRACER_IN

  CHARACTER (LEN=*), INTENT(IN)  :: FILENAME, CATEGORY_IN 

  REAL*8,            INTENT(IN)  :: TAU0_IN
  REAL*4,            INTENT(OUT) :: ARRAY(IX, JX, LX)      

  !**** Local variables
  INTEGER, PARAMETER :: IUNIT = 65
  LOGICAL            :: FOUND 

  INTEGER            :: I,  J,  L,  N,  IOS,M
  INTEGER            :: I1, I2, J1, J2, L1, L2

  REAL*4             :: TEMPARRAY(IIPAR,JJPAR,LLPAR)

  !**** For binary punch file, version 2.0
  INTEGER            :: NTRACER,   NSKIP
  INTEGER            :: HALFPOLAR, CENTER180
  INTEGER            :: NI,        NJ,        NL
  INTEGER            :: IFIRST,    JFIRST,    LFIRST

  REAL*4             :: LONRES,    LATRES

  REAL*8             :: ZTAU0,     ZTAU1

  CHARACTER (LEN=20) :: MODELNAME
  CHARACTER (LEN=40) :: CATEGORY
  CHARACTER (LEN=40) :: UNIT     
  CHARACTER (LEN=40) :: RESERVED
  CHARACTER (LEN=40) :: FTI
  CHARACTER (LEN=80) :: TITLE 

  !
  !*****************************************************************************
  !  READ_BPCH2 begins here!
  !  
  !  Initialize some variables
  !*****************************************************************************
  !
  FOUND            = .FALSE.
  ARRAY(:,:,:)     = 0e0
  TEMPARRAY(:,:,:) = 0e0
  !
  !*****************************************************************************
  !  Open restart file and read header.
  !  Do some error checking to make sure the file is the right format.
  !*****************************************************************************
  !
  
  OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD',&
        FORM='UNFORMATTED',IOSTAT=IOS,CONVERT='BIG_ENDIAN')

  IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:1' )

  ! Read file type identifier
  READ( IUNIT, IOSTAT=IOS ) FTI

  IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:2' )

  IF ( TRIM( FTI ) /= 'CTM bin 02' ) THEN
     PRINT*, 'Input file is not in binary file format v. 2.0!'
     PRINT*, 'STOP in read_bpch2.f'
     STOP
  ENDIF

  ! Read top title
  READ( IUNIT, IOSTAT=IOS ) TITLE

  IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:3' )
  !
  !*****************************************************************************
  !  Read data from the binary punch file 
  !
  !  NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
  !*****************************************************************************
  !
  DO
     READ( IUNIT, IOSTAT=IOS )  MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

     IF ( IOS < 0 ) EXIT
     IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:4' )

     READ( IUNIT, IOSTAT=IOS ) &
             CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,&
             NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,&
             NSKIP


     IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:5' )

     READ( IUNIT, IOSTAT=IOS ) &
             ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

     IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:6' )

     !      print*,  CATEGORY
     !      print*, NTRACER
     !      print*, ZTAU0

     ! Test for a match
     IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and. &
              TRACER_IN           == NTRACER          .and.&
                TAU0_IN             == ZTAU0 ) THEN
        FOUND = .TRUE.
        EXIT
     ENDIF

  ENDDO


  !
  !*****************************************************************************
  !  We have found a match!  Copy TEMPARRAY to ARRAY, taking into account
  !  the starting positions (IFIRST, JFIRST, LFIRST) of the data block.
  !*****************************************************************************
  !
  IF ( FOUND ) THEN 
     I1 = IFIRST
     J1 = JFIRST
     L1 = LFIRST

     I2 = NI + IFIRST - 1
     J2 = NJ + JFIRST - 1
     L2 = NL + LFIRST - 1

     !print*, i1,i2
     !print*, j1,j2
     !print*, L1,L2

     ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )

     !print*, maxval(ARRAY(i1:i2/2,:,:))
     !print*, maxval(array(i2/2+1:i2,:,:))


     !         WRITE( 6, 100 ) ZTAU0, NTRACER
     ! 100     FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2, 
     !     &           ' and tracer # ', i6 )
     !         print*,'NI, NJ, NL, IFIRST, JFIRST, LFIRST = ',
     !     &          NI, NJ, NL, IFIRST, JFIRST, LFIRST 

  ELSE
     WRITE( 6, 110 ) TRIM( FILENAME )
110  FORMAT( 'No matches found for file ', a )
     STOP
  ENDIF
  !
  !*****************************************************************************
  !  Close file and return to calling program  
  !*****************************************************************************
  !
  CLOSE( IUNIT )

END SUBROUTINE READ_BPCH2




