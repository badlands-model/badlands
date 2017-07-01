!     HotConcat
!     Author: Ben Payment 07.17.2006
!
!     Program take hotfiles created from SWAN parallel runs and compiles them
!     into a single file.
!
!     Read command line arguments (only in case of Fortran 2003)
!     HotConcat [-v] [-s] [-h halosize] <basefile>
!     <basefile> refers to the base file name used for hotfiles
!     [-v]           verbose mode: program reports non-error information
!     [-h halosize]  optional argument which changes the overlap halo from default value
!     [-s]           stomp over existing basefile
!     [help] [-help] display above information
!
!     In case of non-Fortran 2003, read arguments via arg_nml namelist from
!     hcat.nml file.
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2011  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!


      program HotConcat

      implicit none
      character(256) RLINE
      character*128 arg_str, basefile, tmp
      integer ios,i,j,k,l,numfiles, sum, halo,haloindex, mxc, myc
      integer prevsource,freq,cdir, version
      logical ydivide, verbose, haloflag, exists, stomp, nonstat, kspher
      real, allocatable :: x(:,:),y(:,:)
      integer, allocatable :: locations(:),index(:), source(:)
      integer, allocatable :: ownlocat(:)
      namelist /arg_nml/ basefile,halo,verbose,stomp

 33   FORMAT('-',I3.3,A)
 35   FORMAT(A)
 36   FORMAT(F14.4,F14.4)
 37   FORMAT(F12.6,F12.6)
 38   FORMAT(I8,2I6,A)

!.....Set command-line argument defaults
      basefile='NULL'
      halo=3
      verbose = .FALSE.
      stomp = .FALSE.

!.....Aquire arguments via arg_nml namelist read from hcat.nml file (non-Fortran 2003)
      OPEN(10,FILE='hcat.nml',STATUS='OLD')
      READ(10,NML=arg_nml)
      CLOSE(10)

!.....Process command-line arguments (Fortran 2003)
!     haloflag = .FALSE.
!     IF(COMMAND_ARGUMENT_COUNT()<1 .OR.
!    &     COMMAND_ARGUMENT_COUNT()>4) THEN
!        PRINT *, 'HotConcat [-v] [-h halosize] <basefile>'
!       RETURN
!     ELSE
!        DO i=1,COMMAND_ARGUMENT_COUNT()
!           CALL GET_COMMAND_ARGUMENT(NUMBER=i,VALUE=arg_str)
!           IF (arg_str == 'help' .OR. arg_str == '-help') THEN
!              PRINT *,'HotConcat [-v] [-s] [-h halosize] <basefile>'
!              PRINT *,'[-v]           verbose mode: program reports ',
!    &                 'non-error information'
!              PRINT *,'[-h halosize]  optional argument which changes',
!    &                 ' the overlap halo from defualt value of 3'
!              PRINT *,'[-s]           stomp over existing basefile'
!              PRINT *,'[help] [-help] display above information'
!              STOP
!           ELSE IF (haloflag) THEN
!              READ (arg_str,'(I10)'), halo
!              haloflag = .FALSE.
!           ELSE IF (arg_str == '-h') THEN
!              haloflag = .TRUE.
!           ELSE IF (arg_str == '-v') THEN
!              verbose = .TRUE.
!           ELSE IF (arg_str == '-s') THEN
!              stomp = .TRUE.
!           ELSE IF (arg_str == '-sv' .OR. arg_str == '-vs') THEN
!              stomp = .TRUE. ; verbose = .TRUE.
!           ELSE IF (arg_str == '-sh' .OR. arg_str == '-hs') THEN
!              stomp = .TRUE. ; haloflag = .TRUE.
!           ELSE IF (arg_str == '-vh' .OR. arg_str == '-hv') THEN
!              verbose = .TRUE. ; haloflag = .TRUE.
!           ELSE IF (arg_str == '-svh' .OR. arg_str == '-shv' .OR.
!    &               arg_str == '-vsh' .OR. arg_str == '-vhs' .OR.
!    &               arg_str == '-hsv' .OR. arg_str == '-hvs') THEN
!              stomp = .TRUE. ; verbose = .TRUE. ;  haloflag = .TRUE.
!           ELSE
!              basefile=arg_str
!           ENDIF
!        ENDDO
!     ENDIF

!.....Error trap
      IF(basefile=='NULL') THEN
         PRINT *,'Invalid basefile or basefile argument missing'
         STOP
      ENDIF

!.....Report halo size
      IF(verbose) PRINT *,'Halosize is',halo

!.....Open basefile
      IF (stomp) THEN
         OPEN(unit=10, file=TRIM(basefile), iostat=ios)
      ELSE
         OPEN(unit=10, file=TRIM(basefile), iostat=ios, STATUS='NEW')
      ENDIF
      IF ( ios /= 0 ) then
         PRINT*,'Could not open ',trim(basefile),' or it already exists'
         STOP
      ENDIF

!.....Identify subdomain files
      DO i=1,1000
         tmp = trim(basefile)
         WRITE(tmp(LEN_TRIM(tmp)+1:LEN_TRIM(tmp)+4),33) i
         INQUIRE(FILE=tmp,EXIST=exists)
         IF(exists) THEN
            numfiles = i
         ELSE
            EXIT
         ENDIF
      ENDDO

!.....Error trap
      IF(numfiles<2) THEN
         PRINT *,'Unable to find more than 1 file'
         STOP
      ENDIF

!.....Report number of files
      IF(verbose) PRINT *,'Trying to open',numfiles,'files'

!.....Open subdomain files
      DO i=1,numfiles
         tmp = trim(basefile)
         WRITE(tmp(LEN_TRIM(tmp)+1:LEN_TRIM(tmp)+4),33) i
         OPEN(unit=10+i, file=tmp, status='old', iostat=ios)
         IF ( ios /= 0 ) THEN
            PRINT *,'Unable to open',TRIM(tmp)
         ELSE
            IF(verbose) PRINT *,TRIM(tmp),' opened succesfully'
         ENDIF
      ENDDO

!.....Specify the number of location for each hotfile including
!.....the one being created at location(0)
      nonstat = .FALSE. !nonstationary unless time is found
      allocate(locations(0:numfiles))
      DO i=1,numfiles
         READ (10+i,35) RLINE
         IF(i==1) WRITE(10,35) TRIM(RLINE)
         READ(RLINE,'(A4,I5)') tmp,version
         IF(tmp /= 'SWAN') THEN
            tmp = trim(basefile)
            WRITE(tmp(LEN_TRIM(tmp)+1:LEN_TRIM(tmp)+4),33) i
            PRINT *,trim(tmp),
     &           'is not a proper hotfile missing header SWAN'
            STOP
         ENDIF
         IF (version /= 1) THEN
            tmp = trim(basefile)
            WRITE(tmp(LEN_TRIM(tmp)+1:LEN_TRIM(tmp)+4),33) i
            PRINT *,trim(tmp),'is verion',version,
     &                'only version 1 is supported'
            STOP
         ENDIF
 110     READ (10+i,35) RLINE
         IF(i==1) WRITE(10,35) TRIM(RLINE)
         IF (RLINE(1:4)=='TIME') nonstat=.TRUE.
         IF (RLINE(1:6).NE.'LONLAT' .AND. RLINE(1:4).NE.'LOCA') GOTO 110
         IF (RLINE(1:6).EQ.'LONLAT') THEN
            kspher = .TRUE.
         ELSE IF (RLINE(1:4).EQ.'LOCA') THEN
            kspher = .FALSE.
         ENDIF
         READ (10+i,38) locations(i),mxc,myc,RLINE
         IF(verbose) PRINT *,'File',i,'has',locations(i),'locations'
      ENDDO

!.....Find sum of location length
      sum=0
      DO i=1,numfiles
         sum=sum+locations(i)
      ENDDO

!.....Allocate 2D array for x,y data
      allocate(x(0:numfiles,1:sum))
      allocate(y(0:numfiles,1:sum))

!.....Allocate array to refer to which hotfile a grid point info is stored
      allocate(source(1:sum))

!.....Read in locations
      IF (.not.kspher) THEN
         DO i=1,numfiles
            DO j=1,locations(i)
               READ (10+i,36) x(i,j),y(i,j)
            ENDDO
         ENDDO
      ELSE
         DO i=1,numfiles
            DO j=1,locations(i)
               READ (10+i,37) x(i,j),y(i,j)
            ENDDO
         ENDDO
      ENDIF

!.....Determine how data is divided for parallel run either on x or y
      IF(x(1,1)==x(2,1)) THEN
         ydivide=.TRUE.
         IF(verbose) PRINT *,'Processes divided along y'
      ELSE
         ydivide=.FALSE.
         IF(verbose) PRINT *,'Processes divided along x'
      ENDIF

!.....Index array keep track of the index number of the grid point being processes
      allocate(index(0:numfiles))
      DO i=0,numfiles
         index(i)=1
      ENDDO

!.....Sort grid points
      DO WHILE (index(numfiles)<=locations(numfiles))
!........IF divided along Y
         IF(ydivide) THEN
            DO i=1,numfiles
               IF(i==numfiles) THEN
                  DO WHILE(x(i,index(i))<=x(1,index(1)-1))
                        x(0,index(0))=x(i,index(i))
                        y(0,index(0))=y(i,index(i))
                        source(index(0))=i
                        index(i)=index(i)+1
                        index(0)=index(0)+1
                        IF(index(i)>locations(i)) EXIT
                  ENDDO
               ELSEIF (index(i)<=locations(i)) THEN
                  haloindex=0
                  DO WHILE(x(i,index(i))<=x(i+1,index(i+1)) .AND.
     &                       haloindex<halo)
                     x(0,index(0))=x(i,index(i))
                     y(0,index(0))=y(i,index(i))
                     source(index(0))=i
                     IF( y(i,index(i)) == y(i+1,index(i+1)) ) THEN
                        haloindex=haloindex+1
                        index(i+1)=index(i+1)+1
                     ENDIF
                     index(i)=index(i)+1
                     index(0)=index(0)+1
                  ENDDO
                  index(i)=index(i)+halo
               ELSE
                  PRINT *,'Error'
                  STOP
               ENDIF
            ENDDO
!........IF divided along X
         ELSE
            DO i=1,numfiles
               IF(i==numfiles) THEN
                  DO WHILE(index(i)<=locations(i))
                     x(0,index(0))=x(i,index(i))
                     y(0,index(0))=y(i,index(i))
                     source(index(0))=i
                     index(i)=index(i)+1
                     index(0)=index(0)+1
                  ENDDO
               ELSE
                  haloindex=0
                  DO WHILE(x(i,index(i))<=x(i+1,index(i+1)) .AND.
     &                 haloindex<halo)
                     x(0,index(0))=x(i,index(i))
                     y(0,index(0))=y(i,index(i))
                     source(index(0))=i
                     IF( x(i,index(i)) == x(i+1,index(i+1)) .AND.
     &                   y(i,index(i)) == y(i+1,index(i+1)) ) THEN
                        IF (y(i+1,index(i+1))>y(i+1,index(i+1)+1)) THEN
                           haloindex=haloindex+1
                        ENDIF
                        index(i+1)=index(i+1)+1
                     ENDIF
                     index(i)=index(i)+1
                     index(0)=index(0)+1
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!.....Index shift, set total number of locations
      DO i=0,numfiles
         index(i)=index(i)-1
      ENDDO
      locations(0)=index(0)

!.....Specify the number of own location for each hotfile
      allocate(ownlocat(0:numfiles))
      ownlocat=0
      DO j=1,numfiles
         index(0)=1
         DO WHILE (index(0)<=locations(0))
            IF (source(index(0))==j) ownlocat(j)=ownlocat(j)+1
            index(0)=index(0)+1
         END DO
      END DO

      IF (ydivide) THEN
         myc = locations(0)/mxc
      ELSE
         mxc = locations(0)/myc
      ENDIF

!.....Write total number of locations to basefile
      WRITE (10,38) locations(0),mxc,myc,TRIM(RLINE)

!.....Write grid points to basefile
      IF (.not.kspher) THEN
         DO i=1,locations(0)
            WRITE (10,36) x(0,i),y(0,i)
         ENDDO
      ELSE
         DO i=1,locations(0)
            WRITE (10,37) x(0,i),y(0,i)
         ENDDO
      ENDIF

!.....Process FREQ
      DO i=1,numfiles
         READ (10+i,35) RLINE
         IF (i==1) WRITE(10,35) TRIM(RLINE)
         READ (10+i,35) RLINE
         READ (RLINE,*) freq
         IF (i==1) WRITE(10,35) TRIM(RLINE)
         IF(verbose) PRINT *,'File',i,'has',freq,'frequencies'
         DO j=1,freq
            READ (10+i,35) RLINE
            IF (i==1) WRITE(10,35) TRIM(RLINE)
         ENDDO
      ENDDO

!.....Process DIR
      DO i=1,numfiles
         READ (10+i,35) RLINE
         IF (i==1) WRITE(10,35) TRIM(RLINE)
         READ (10+i,35) RLINE
         READ (RLINE,*) cdir
         IF (i==1) WRITE(10,35) TRIM(RLINE)
         IF(verbose) PRINT *,'File',i,'has',cdir,'directions'
         DO j=1,cdir
            READ (10+i,35) RLINE
            IF (i==1) WRITE(10,35) TRIM(RLINE)
         ENDDO
      ENDDO

!.....Process QUANT (must be single quantity)
      DO i=1,numfiles
         DO j=1,5
            READ (10+i,35) RLINE
            IF (i==numfiles) WRITE(10,35) TRIM(RLINE)
         ENDDO
      ENDDO

!.....Process date & time -- in case of nonstationary run
      IF (nonstat) THEN
         DO i=1,numfiles
            READ (10+i,35) RLINE
            IF (i==numfiles) WRITE(10,35) TRIM(RLINE)
         ENDDO
      ENDIF

!.....Reset index
      DO i=0,numfiles
         index(i)=1
      ENDDO

!.....Change in source represent hotfile boundary and the need to process halo points
      prevsource=1

!.....Process action densities
      DO WHILE (index(0)<=locations(0))
         DO j=1,numfiles
!...........IF divided along Y
            IF(ydivide) THEN
               IF(source(index(0))==j) THEN
!.................IF source equals 1 there are no halo points at the beginning
                  IF(j .NE. 1) THEN
                     DO k=1,halo
                        READ (10+j,35,IOSTAT=ios) RLINE
                        IF (RLINE(1:6).EQ.'FACTOR' .and. ios==0) THEN
                           DO l=1,freq+1
                              READ (10+j,35,IOSTAT=ios) RLINE
                           ENDDO
                           index(j)=index(j)+1
                        ELSE IF(RLINE(1:4).EQ.'ZERO' .and. ios==0) THEN
                           index(j)=index(j)+1
                        ELSE IF(RLINE(1:6).EQ.'NODATA' .and. ios==0)THEN
                           index(j)=index(j)+1
                        ELSE
                           PRINT *,'Unable to process'
                           PRINT *,RLINE
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
!.................Write out grid point data as long as the source doesn't change
                  DO WHILE(source(index(0))==j)
                     READ (10+j,35,IOSTAT=ios) RLINE
                     IF (RLINE(1:6).EQ.'FACTOR' .and. ios==0) THEN
                        WRITE(10,35) TRIM(RLINE)
                        DO l=1,freq+1
                           READ (10+j,35,IOSTAT=ios) RLINE
                           WRITE(10,35) TRIM(RLINE)
                        ENDDO
                        index(j)=index(j)+1
                        index(0)=index(0)+1
                     ELSE IF(RLINE(1:4).EQ.'ZERO' .and. ios==0) THEN
                        WRITE(10,35) TRIM(RLINE)
                        index(j)=index(j)+1
                        index(0)=index(0)+1
                     ELSE IF(RLINE(1:6).EQ.'NODATA' .and. ios==0) THEN
                        WRITE(10,35) TRIM(RLINE)
                        index(j)=index(j)+1
                        index(0)=index(0)+1
                     ELSE
                        PRINT *,'Unable to process'
                        PRINT *,RLINE
                        EXIT
                     ENDIF
                     prevsource=j
                  ENDDO
!.................IF source equals last source there are no halo points at the end
                  IF(j .NE. numfiles) THEN
                     DO k=1,halo
                        READ (10+j,35,IOSTAT=ios) RLINE
                        IF (RLINE(1:6).EQ.'FACTOR' .and. ios==0) THEN
                           DO l=1,freq+1
                              READ (10+j,35,IOSTAT=ios) RLINE
                           ENDDO
                           index(j)=index(j)+1
                        ELSE IF(RLINE(1:4).EQ.'ZERO' .and. ios==0) THEN
                           index(j)=index(j)+1
                        ELSE IF(RLINE(1:6).EQ.'NODATA' .and. ios==0)THEN
                           index(j)=index(j)+1
                        ELSE
                           PRINT *,'Unable to process'
                           PRINT *,RLINE
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
!...........IF divided along X
            ELSE
               IF(source(index(0))==j) THEN
!.................IF source equals 1 there are no halo points at the beginning
                  IF(j .NE. 1) THEN
                     DO WHILE (index(j-1)<=ownlocat(j-1))
                        READ (10+j,35,IOSTAT=ios) RLINE
                        IF (RLINE(1:6).EQ.'FACTOR' .and. ios==0)THEN
                           DO l=1,freq+1
                              READ (10+j,35,IOSTAT=ios) RLINE
                           ENDDO
                           index(j-1)=index(j-1)+1
                           index(j)=index(j)+1
                        ELSEIF(RLINE(1:4).EQ.'ZERO'.and.ios==0) THEN
                           index(j-1)=index(j-1)+1
                           index(j)=index(j)+1
                        ELSEIF(RLINE(1:6).EQ.'NODATA'.and.ios==0)THEN
                           index(j-1)=index(j-1)+1
                           index(j)=index(j)+1
                        ELSE IF (ios/=0) THEN
                           PRINT *,'End of file error'
                           PRINT *,'ios:',ios
                           STOP
                        ELSE
                           PRINT *,'Unable to process'
                           PRINT *,'RLINE: ',RLINE
                           STOP
                        ENDIF
                     ENDDO
                  ENDIF
!.................Write out grid point data as long as the source doesn't change
                  DO WHILE(source(index(0))==j)
                     READ (10+j,35,IOSTAT=ios) RLINE
                     IF (RLINE(1:6).EQ.'FACTOR' .and. ios==0) THEN
                        WRITE(10,35) TRIM(RLINE)
                        DO l=1,freq+1
                           READ (10+j,35,IOSTAT=ios) RLINE
                           WRITE(10,35) TRIM(RLINE)
                        ENDDO
                        index(j)=index(j)+1
                        index(0)=index(0)+1
                     ELSE IF(RLINE(1:4).EQ.'ZERO' .and. ios==0) THEN
                        WRITE(10,35) TRIM(RLINE)
                        index(j)=index(j)+1
                        index(0)=index(0)+1
                     ELSE IF(RLINE(1:6).EQ.'NODATA' .and. ios==0) THEN
                        WRITE(10,35) TRIM(RLINE)
                        index(j)=index(j)+1
                        index(0)=index(0)+1
                     ELSE IF (ios/=0) THEN
                        PRINT *,'End of file error'
                        PRINT *,'ios:',ios
                        STOP
                     ELSE
                        PRINT *,'Unable to process'
                        PRINT *,'RLINE: ',RLINE
                        STOP
                     ENDIF
                     prevsource=j
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!.....Close all files
      CLOSE ( unit = 10 );
      DO i=1,numfiles
         CLOSE (unit=10+i)
      ENDDO

      END PROGRAM HotConcat

