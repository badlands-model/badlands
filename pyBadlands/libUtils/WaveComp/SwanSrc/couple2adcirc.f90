      MODULE Couple2Adcirc

      INTEGER   IUNIT
      INTEGER   IOSTAT, IT0, SAVITE, ILEN
      INTEGER   INERR
      INTEGER   IERR
      INTEGER   ISTAT, IF1, IL1
      CHARACTER PTYPE, PNAME *8, COMPUT *4
      CHARACTER*20 CHARS(1)
      CHARACTER*20 MSGSTR
      LOGICAL   LOPEN

      INTEGER, ALLOCATABLE :: CROSS(:)
      INTEGER, ALLOCATABLE :: BGRIDP(:)
      REAL   , ALLOCATABLE :: BSPECS(:,:,:,:)
      REAL   , ALLOCATABLE :: AC1(:,:,:), COMPDA(:,:)

      REAL, ALLOCATABLE    :: BLKND(:), BLKNDC(:), OURQT(:)

      CONTAINS



      SUBROUTINE MakeBoundariesReflective(ivert, ac2 )

      USE M_GENARR, ONLY        : spcdir
      USE SwanGridData, ONLY    : nverts
      USE SwanGridobjects, ONLY : CELLID,               &
                                  celltype,             &
                                  CELLV1,CELLV2,CELLV3, &
                                  FACEID,               &
                                  facetype,             &
                                  FACEV1,FACEV2,        &
                                  FMARKER,              &
                                  gridobject,           &
                                  verttype,             &
                                  VERTX,VERTY
      USE SWCOMM3, ONLY         : DDIR,                 &
                                  MDC,MSC
      IMPLICIT NONE

      INTRINSIC              :: ATAN
      INTRINSIC              :: COS
      INTRINSIC              :: INT
      INTRINSIC              :: MOD
      INTRINSIC              :: REAL
      INTRINSIC              :: SIN
      INTRINSIC              :: SQRT

      INTEGER                :: icell
      INTEGER                :: id
      INTEGER                :: iface
      INTEGER                :: IncD1
      INTEGER                :: IncD2
      INTEGER                :: is
      INTEGER,INTENT(IN)     :: ivert
      INTEGER                :: jc
      INTEGER                :: jf
      INTEGER                :: NumBdySegs
      INTEGER                :: NumFaces
      INTEGER                :: V1
      INTEGER                :: V2
      INTEGER                :: V3

      LOGICAL                :: IntoDomain

      REAL,INTENT(INOUT)     :: ac2(MDC,MSC,nverts)
      REAL                   :: ac2temp(MDC,MSC)
      REAL                   :: AvgDX
      REAL                   :: BdyDir(2)
      REAL                   :: BdyDirTemp
      REAL                   :: BdyDirToUse
      REAL                   :: Dir
      REAL                   :: Dist(2)
      REAL                   :: IncDir
      REAL                   :: IncCounter
      REAL                   :: Pi = 3.141592654
      REAL                   :: SubArea1
      REAL                   :: SubArea2
      REAL                   :: SubArea3
      REAL                   :: TotalArea
      REAL                   :: W1
      REAL                   :: W2
      REAL                   :: X1
      REAL                   :: X2
      REAL                   :: X3
      REAL                   :: XV
      REAL                   :: Y1
      REAL                   :: Y2
      REAL                   :: Y3
      REAL                   :: YV

      TYPE(celltype),POINTER :: cell(:)
      TYPE(facetype),POINTER :: face(:)
      TYPE(verttype),POINTER :: vert(:)

!... Initialize ac2temp.
      DO id=1,MDC
         DO is=1,MSC
            ac2temp(id,is) = 0.
         ENDDO
      ENDDO

!... Point to data.
      cell => gridobject%cell_grid
      face => gridobject%face_grid
      vert => gridobject%vert_grid

!... Find a length scale by taking the average dx between this
!... vertex and its neighbors.
      AvgDX = 0.
      NumFaces = 0
      DO jc=1,vert(ivert)%noc
         icell = vert(ivert)%cell(jc)%atti(CELLID)
         DO jf=1,cell(icell)%nof
            V1 = cell(icell)%face(jf)%atti(FACEV1)
            V2 = cell(icell)%face(jf)%atti(FACEV2)
            IF( V1==ivert .OR. V2==ivert )THEN
               X1 = vert(V1)%attr(VERTX)
               X2 = vert(V2)%attr(VERTX)
               Y1 = vert(V1)%attr(VERTY)
               Y2 = vert(V2)%attr(VERTY)
               AvgDX = AvgDX + SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))
               NumFaces = NumFaces + 1
            ENDIF
         ENDDO
      ENDDO
      AvgDX = AvgDX / REAL(NumFaces)

!... Loop over the directions.
      DO id=1,MDC
         Dir = spcdir(id,1)

!... Determine if this direction points out of the domain.

!... Set up a fictional vertex along the direction.
         XV = vert(ivert)%attr(VERTX) + COS(Dir) * ( 0.1 * AvgDX )
         YV = vert(ivert)%attr(VERTY) + SIN(Dir) * ( 0.1 * AvgDX )

!... Loop over the adjacent cells.  Use triangle areas
!... to determine if the fictional vertex is located
!... inside an adjacent cell.
         IntoDomain = .FALSE.
         DO jc=1,vert(ivert)%noc

            icell = vert(ivert)%cell(jc)%atti(CELLID)
            V1 = cell(icell)%atti(CELLV1)
            V2 = cell(icell)%atti(CELLV2)
            V3 = cell(icell)%atti(CELLV3)

            X1 = XV
            X2 = vert(V2)%attr(VERTX)
            X3 = vert(V3)%attr(VERTX)
            Y1 = YV
            Y2 = vert(V2)%attr(VERTY)
            Y3 = vert(V3)%attr(VERTY)
            SubArea1 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            X1 = vert(V1)%attr(VERTX)
            X2 = XV
            X3 = vert(V3)%attr(VERTX)
            Y1 = vert(V1)%attr(VERTY)
            Y2 = YV
            Y3 = vert(V3)%attr(VERTY)
            SubArea2 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            X1 = vert(V1)%attr(VERTX)
            X2 = vert(V2)%attr(VERTX)
            X3 = XV
            Y1 = vert(V1)%attr(VERTY)
            Y2 = vert(V2)%attr(VERTY)
            Y3 = YV
            SubArea3 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            X1 = vert(V1)%attr(VERTX)
            X2 = vert(V2)%attr(VERTX)
            X3 = vert(V3)%attr(VERTX)
            Y1 = vert(V1)%attr(VERTY)
            Y2 = vert(V2)%attr(VERTY)
            Y3 = vert(V3)%attr(VERTY)
            TotalArea = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            IF( (SubArea1+SubArea2+SubArea3).LE.(1.001*TotalArea) )THEN

               IntoDomain = .TRUE.

            ENDIF

         ENDDO

!... If the fictional vertex does not lie inside the domain,
!... then we know that the current direction is pointing outward.
!... Reflect its wave energy back into the domain.
         IF( .NOT. IntoDomain )THEN

!... Find a representative direction for the boundary by taking
!... an average of the boundary segment(s) connected to this vertex.

!... Loop over adjacent faces to find boundary segments.
            NumBdySegs = 0
            DO jc=1,vert(ivert)%noc
               icell = vert(ivert)%cell(jc)%atti(CELLID)

               DO jf=1,cell(icell)%nof
                  iface = cell(icell)%face(jf)%atti(FACEID)
                  V1 = face(iface)%atti(FACEV1)
                  V2 = face(iface)%atti(FACEV2)

!... If this face is a boundary segment, then compute its direction.
                  IF( face(iface)%atti(FMARKER)==1 .AND. &
                     (V1==ivert .OR. V2==ivert) )THEN

                     X1 = vert(ivert)%attr(VERTX)
                     Y1 = vert(ivert)%attr(VERTY)
                     IF( V1==ivert )THEN
                        X2 = vert(V2)%attr(VERTX)
                        Y2 = vert(V2)%attr(VERTY)
                     ELSE
                        X2 = vert(V1)%attr(VERTX)
                        Y2 = vert(V1)%attr(VERTY)
                     ENDIF

                     NumBdySegs = NumBdySegs + 1
                     BdyDirTemp = ATAN( (Y2-Y1) / (X2-X1) )
                     IF( (Y2-Y1).GE.0. .AND. (X2-X1).GE.0. )THEN
                        BdyDirTemp = BdyDirTemp
                     ELSEIF( (Y2-Y1).GE.0. .AND. (X2-X1).LT.0. )THEN
                        BdyDirTemp = Pi + BdyDirTemp
                     ELSEIF( (Y2-Y1).LT.0. .AND. (X2-X1).LT.0. )THEN
                        BdyDirTemp = Pi + BdyDirTemp
                     ELSE
                        BdyDirTemp = (2.*Pi) + BdyDirTemp
                     ENDIF
                     BdyDir(NumBdySegs) = BdyDirTemp

                  ENDIF

               ENDDO

            ENDDO

!... Determine the boundary direction across which to reflect
!... the wave energy.

!... Determine the distances between the current direction
!... and the directions of the two boundary segments.
            Dist(1) = ABS(BdyDir(1)-Dir)
            IF( Dist(1).GT.Pi )THEN
               Dist(1) = ABS(BdyDir(1)-(Dir-(2.*Pi)))
            ENDIF
            Dist(2) = ABS(BdyDir(2)-Dir)
            IF( Dist(2).GT.Pi )THEN
               Dist(2) = ABS(BdyDir(2)-(Dir-(2.*Pi)))
            ENDIF

!... If the current direction is more than 90 degrees
!... from both of the boundary segments, then reflect
!... its wave energy across the average direction of the
!... boundary segments.
            IF( Dist(1).GT.(0.5*Pi) .AND. &
                Dist(2).GT.(0.5*Pi) )THEN

               BdyDirToUse = ( (BdyDir(1)-Pi) + &
                               BdyDir(2) ) / 2.0
               IF( BdyDirToUse.LT.0. )THEN
                  BdyDirToUse = BdyDirToUse + (2.*Pi)
               ENDIF

!... If the current direction is an equal distance
!... from both of the boundary segments, then reflect
!... its wave energy across the average direction of the
!... boundary segments.
            ELSEIF( Dist(1).EQ.Dist(2) )THEN

               BdyDirToUse = ( (BdyDir(1)-Pi) + &
                               BdyDir(2) ) / 2.0
               IF( BdyDirToUse.LT.0. )THEN
                  BdyDirToUse = BdyDirToUse + (2.*Pi)
               ENDIF

!... Otherwise, reflect the wave energy across
!... the closest boundary segment.
            ELSE

               IF( Dist(1).LT.Dist(2) )THEN
                  BdyDirToUse = BdyDir(1)
               ELSE
                  BdyDirToUse = BdyDir(2)
               ENDIF

            ENDIF

!... Reflect our current direction across the boundary direction.
            IncDir = 2. * BdyDirToUse - Dir
            IF( IncDir<0. )THEN
               IncDir = IncDir + (2.*Pi)
            ENDIF

!... Determine counter for which direction is IncDir.
            IncCounter = MOD( IncDir - spcdir(1,1), (2.*Pi) ) / DDIR
            IF( IncCounter<0. )THEN
               IncCounter = IncCounter + REAL(MDC)
            ENDIF

!... Determine the directions in our spectrum to which we can reflect
!... wave energy, and how much energy should go to each.
            IncD1 = 1 + INT(IncCounter)
            IncD2 = IncD1 + 1
            IF( IncD2>MDC )THEN
               IncD2 = IncD2 - MDC
            ENDIF
            W2 = IncCounter + 1. - REAL(IncD1)
            W1 = 1. - W2

!... It is sometimes possible for the incident direction to fall between
!... two directions that are not inside the active domain.  In this case,
!... one direction will be at the edge of the active domain (and near
!... the boundary), and the other direction will be just outside
!... of the domain.  We don't want to reflect energy onto that latter direction.
!... Adjust the weights to take care of this case.
            IF( IncD1==id .OR. IncD2==id )THEN
               IF( IncD1==id )THEN
                  W1 = 0.
                  W2 = 1.
               ELSE
                  W1 = 1.
                  W2 = 0.
               ENDIF
            ENDIF

!... Loop over spectrum, reflecting the wave energy as appropriate.
            DO is=1,MSC

!... Reflect components of wave energy into the domain.
               ac2temp(IncD1,is) = W1 * ac2(id,is,ivert) + ac2temp(IncD1,is)
               ac2temp(IncD2,is) = W2 * ac2(id,is,ivert) + ac2temp(IncD2,is)

            ENDDO

         ENDIF

      ENDDO

!... We must be very careful about how we add the reflected energy
!... back into the computational domain.  Our method must depend
!... on whether the direction: (1) points out of the domain,
!... (2) points into the domain through the boundary, or
!... (3) points into the domain from inside the domain.

!... First, determine if the current direction points
!... out of the domain.  Use the same method as above.
      DO id=1,MDC
         Dir = spcdir(id,1)

!... Set up a fictional vertex along the direction.
         XV = vert(ivert)%attr(VERTX) + COS(Dir) * ( 0.1 * AvgDX )
         YV = vert(ivert)%attr(VERTY) + SIN(Dir) * ( 0.1 * AvgDX )

!... Loop over the adjacent cells.  Use triangle areas
!... to determine if the fictional vertex is located
!... inside an adjacent cell.
         IntoDomain = .FALSE.
         DO jc=1,vert(ivert)%noc

            icell = vert(ivert)%cell(jc)%atti(CELLID)
            V1 = cell(icell)%atti(CELLV1)
            V2 = cell(icell)%atti(CELLV2)
            V3 = cell(icell)%atti(CELLV3)

            X1 = XV
            X2 = vert(V2)%attr(VERTX)
            X3 = vert(V3)%attr(VERTX)
            Y1 = YV
            Y2 = vert(V2)%attr(VERTY)
            Y3 = vert(V3)%attr(VERTY)
            SubArea1 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            X1 = vert(V1)%attr(VERTX)
            X2 = XV
            X3 = vert(V3)%attr(VERTX)
            Y1 = vert(V1)%attr(VERTY)
            Y2 = YV
            Y3 = vert(V3)%attr(VERTY)
            SubArea2 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            X1 = vert(V1)%attr(VERTX)
            X2 = vert(V2)%attr(VERTX)
            X3 = XV
            Y1 = vert(V1)%attr(VERTY)
            Y2 = vert(V2)%attr(VERTY)
            Y3 = YV
            SubArea3 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            X1 = vert(V1)%attr(VERTX)
            X2 = vert(V2)%attr(VERTX)
            X3 = vert(V3)%attr(VERTX)
            Y1 = vert(V1)%attr(VERTY)
            Y2 = vert(V2)%attr(VERTY)
            Y3 = vert(V3)%attr(VERTY)
            TotalArea = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

            IF( (SubArea1+SubArea2+SubArea3).LE.(1.001*TotalArea) )THEN

               IntoDomain = .TRUE.

            ENDIF

         ENDDO

!... If the direction points out of the domain, then simply allow
!... its wave energy to remain constant.  Do not change its action
!... densities at all.

         IF( .NOT. IntoDomain )THEN

            CONTINUE

         ELSE

!... If the direction points into the domain, then we must determine
!... if it is coming from inside or outside of the domain.
!... Use the same method as before.
            Dir = spcdir(id,1)
            IF( Dir.LT.Pi )THEN
               Dir = Dir + Pi
            ELSE
               Dir = Dir - Pi
            ENDIF

!... Set up a fictional vertex along the direction.
            XV = vert(ivert)%attr(VERTX) + COS(Dir) * ( 0.1 * AvgDX )
            YV = vert(ivert)%attr(VERTY) + SIN(Dir) * ( 0.1 * AvgDX )

!... Loop over the adjacent cells.  Use triangle areas
!... to determine if the fictional vertex is located
!... inside an adjacent cell.
            IntoDomain = .FALSE.
            DO jc=1,vert(ivert)%noc

               icell = vert(ivert)%cell(jc)%atti(CELLID)
               V1 = cell(icell)%atti(CELLV1)
               V2 = cell(icell)%atti(CELLV2)
               V3 = cell(icell)%atti(CELLV3)

               X1 = XV
               X2 = vert(V2)%attr(VERTX)
               X3 = vert(V3)%attr(VERTX)
               Y1 = YV
               Y2 = vert(V2)%attr(VERTY)
               Y3 = vert(V3)%attr(VERTY)
               SubArea1 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

               X1 = vert(V1)%attr(VERTX)
               X2 = XV
               X3 = vert(V3)%attr(VERTX)
               Y1 = vert(V1)%attr(VERTY)
               Y2 = YV
               Y3 = vert(V3)%attr(VERTY)
               SubArea2 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

               X1 = vert(V1)%attr(VERTX)
               X2 = vert(V2)%attr(VERTX)
               X3 = XV
               Y1 = vert(V1)%attr(VERTY)
               Y2 = vert(V2)%attr(VERTY)
               Y3 = YV
               SubArea3 = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

               X1 = vert(V1)%attr(VERTX)
               X2 = vert(V2)%attr(VERTX)
               X3 = vert(V3)%attr(VERTX)
               Y1 = vert(V1)%attr(VERTY)
               Y2 = vert(V2)%attr(VERTY)
               Y3 = vert(V3)%attr(VERTY)
               TotalArea = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))

               IF( (SubArea1+SubArea2+SubArea3).LE.(1.001*TotalArea) )THEN

                  IntoDomain = .TRUE.

               ENDIF

            ENDDO

!... If the direction is coming from outside of the domain,
!... then it is not receiving any natural wave energy as part
!... of the overall Swan computations.  No wave energy can propagate
!... through the boundary and into the direction.  Thus, we can
!... simply copy the reflected wave energy into these spectral bins.

            IF( .NOT. IntoDomain )THEN

               DO is=1,MSC
                  ac2(id,is,ivert) = ac2temp(id,is)
               ENDDO

!... If the direction is coming from inside the domain,
!... then it will be receiving wave energy in addition to the energy
!... that we are reflecting.  We cannot simply assign the reflected
!... energy into the spectral bins, because we must consider the
!... other wave energy.  However, we cannot simply add the reflected
!... energy, because we only want to reflect the energy once per
!... iteration.

!... So let's only reflect the energy if the wave energy
!... in this direction has just been updated.

            ELSE

               DO is=1,MSC
!                 ac2(id,is,ivert) = ac2(id,is,ivert) + ac2temp(id,is)
                  ac2(id,is,ivert) = ac2temp(id,is)
               ENDDO

            ENDIF

         ENDIF

      ENDDO

      END SUBROUTINE



      END MODULE Couple2Adcirc
