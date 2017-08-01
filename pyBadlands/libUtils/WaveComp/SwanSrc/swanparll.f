!
!     SWAN - routines for distributed-memory approach based on MPI
!
!  Contents of this file
!
!     SWINITMPI
!     SWEXITMPI
!     SWSYNC
!     SWSENDNB
!     SWRECVNB
!     SWBROADC
!     SWGATHER
!     SWREDUCE
!     SWREDUCI
!     SWREDUCR
!     SWSTRIP
!     SWPARTIT
!     SWBLADM
!     SWDECOMP
!     SWEXCHG
!     SWRECVAC
!     SWSENDAC
!     SWCOLLECT
!     SWCOLOUT
!     SWCOLTAB
!     SWCOLSPC
!     SWCOLBLK
!
!****************************************************************
!
      SUBROUTINE SWINITMPI( MY_COMM )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE

      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Join parallel application
!
!  3. Method
!
!     Start MPI and initialize some variables
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IERR, IF1, IF2, IL1, IL2, MY_COMM
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_COMM_RANK    Get rank of processes in MPI communication context
!     MPI_COMM_SIZE    Get number of processes in MPI communication context
!     MPI_INIT         Enroll in MPI
!     MSGERR           Writes error message
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Start MPI and initialize some common variables in module M_PARALL
!
! 13. Source text
!
      LEVERR = 0
      MAXERR = 1
      ITRACE = 0
      WAV_COMM_WORLD = MY_COMM

!     --- enroll in MPI

!      CALL MPI_INIT ( IERR )
!      IF ( IERR.NE.MPI_SUCCESS ) THEN
!         CHARS(1) = INTSTR(IERR)
!         CALL TXPBLA(CHARS(1),IF1,IL1)
!         MSGSTR = 'MPI produces some internal error - '//
!     &            'return code is '//CHARS(1)(IF1:IL1)
!         CALL MSGERR ( 4, MSGSTR )
!         RETURN
!      END IF

!     --- initialize common variables

      INODE = 0
      NPROC = 1

!     --- get node number INODE

      CALL MPI_COMM_RANK ( WAV_COMM_WORLD, INODE, IERR )
      INODE = INODE + 1
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS(1) = INTSTR(IERR)
         CALL TXPBLA(CHARS(1),IF1,IL1)
         CHARS(2) = INTSTR(INODE)
         CALL TXPBLA(CHARS(2),IF2,IL2)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(1)(IF1:IL1)//
     &            ' and node number is '//CHARS(2)(IF2:IL2)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

!     --- determine total number of processes

      CALL MPI_COMM_SIZE ( WAV_COMM_WORLD, NPROC, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS(1) = INTSTR(IERR)
         CALL TXPBLA(CHARS(1),IF1,IL1)
         CHARS(2) = INTSTR(INODE)
         CALL TXPBLA(CHARS(2),IF2,IL2)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(1)(IF1:IL1)//
     &            ' and node number is '//CHARS(2)(IF2:IL2)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

!     --- determine whether this is a parallel run or not

      IF ( NPROC.GT.1 ) THEN
         PARLL = .TRUE.
      ELSE
         PARLL = .FALSE.
      END IF

!     --- am I master?

      IAMMASTER = INODE.EQ.MASTER

!     --- define MPI constants for communication within SWAN

      SWINT  = MPI_INTEGER
      SWREAL = MPI_REAL
      SWMAX  = MPI_MAX
      SWMIN  = MPI_MIN
      SWSUM  = MPI_SUM

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXITMPI
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Exit parallel application
!
!  3. Method
!
!     Wrapper for MPI_FINALIZE
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IERR    :   error value of MPI call
!     PARALMPI:   if true, parallel process is carried out with MPI
!
      INTEGER IERR
      LOGICAL PARALMPI
!
!  8. Subroutines used
!
!     MPI_ABORT        Abort MPI if severe error occurs
!     MPI_BARRIER      Blocks until all nodes have called this routine
!     MPI_INITIALIZED  Indicates whether MPI_Init has been called
!     MPI_FINALIZE     Cleans up the MPI state and exits
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     if MPI has been initialized
!        synchronize nodes
!        if severe error
!           abort MPI
!        else
!           close MPI
!
! 13. Source text
!
      CALL MPI_INITIALIZED ( PARALMPI, IERR )
      IF ( PARALMPI ) THEN

         CALL MPI_BARRIER ( WAV_COMM_WORLD, IERR )

         IF ( LEVERR.GE.4 ) THEN

!        --- in case of a severe error abort all MPI processes

            CALL MPI_ABORT ( WAV_COMM_WORLD, LEVERR, IERR )

         ELSE

!        --- otherwise stop MPI operations on this computer

!            CALL MPI_FINALIZE ( IERR )

         END IF

      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSYNC
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Synchronize nodes
!
!  3. Method
!
!     Wrapper for MPI_BARRIER
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_BARRIER      Blocks until all nodes have called this routine
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Blocks until all nodes have called MPI_BARRIER routine.
!     In this way, all nodes are synchronized
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSYNC')

!     --- blocks until all nodes have called this routine

      CALL MPI_BARRIER ( WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS(1) = INTSTR(IERR)
         CALL TXPBLA(CHARS(1),IF1,IL1)
         CHARS(2) = INTSTR(INODE)
         CALL TXPBLA(CHARS(2),IF2,IL2)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(1)(IF1:IL1)//
     &            ' and node number is '//CHARS(2)(IF2:IL2)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDNB ( IPTR, ILEN, ITYPE, IDEST, ITAG )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Data is sent to a neighbour
!
!  3. Method
!
!     Wrapper for MPI_SEND
!
!  4. Argument variables
!
!     IDEST       rank of the destination process
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, IDEST, ITAG
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_SEND         Immediately sends the data in the active
!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Data is sent to a neighbour with command MPI_SEND
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      CALL MPI_SEND ( IPTR, ILEN, ITYPE, IDEST-1,
     &                ITAG, WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS(1) = INTSTR(IERR)
         CALL TXPBLA(CHARS(1),IF1,IL1)
         CHARS(2) = INTSTR(INODE)
         CALL TXPBLA(CHARS(2),IF2,IL2)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(1)(IF1:IL1)//
     &            ' and node number is '//CHARS(2)(IF2:IL2)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVNB ( IPTR, ILEN, ITYPE, ISOURCE, ITAG )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Data is received from a neighbour
!
!  3. Method
!
!     Wrapper for MPI_RECV
!
!  4. Argument variables
!
!     ILEN        length of array to be received
!     IPTR        pointer to first element of array to be received
!     ISOURCE     rank of the source process
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, ISOURCE, ITAG
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     ISTAT :     MPI status array
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      INTEGER      ISTAT(MPI_STATUS_SIZE)
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_RECV         Immediately receives the data in the active
!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Data is received from a neighbour with command MPI_RECV
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      CALL MPI_RECV ( IPTR, ILEN, ITYPE, ISOURCE-1, ITAG,
     &                WAV_COMM_WORLD, ISTAT, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS(1) = INTSTR(IERR)
         CALL TXPBLA(CHARS(1),IF1,IL1)
         CHARS(2) = INTSTR(INODE)
         CALL TXPBLA(CHARS(2),IF2,IL2)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(1)(IF1:IL1)//
     &            ' and node number is '//CHARS(2)(IF2:IL2)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBROADC ( IPTR, ILEN, ITYPE )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Broadcasts data from the master to all other processes
!
!  3. Method
!
!     Wrapper for MPI_BCAST
!
!  4. Argument variables
!
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_BCAST        Broadcasts a message from the master
!                      to all other processes of the group
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SNEXTI
!     FLFILE
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Broadcasts data from the master to all other nodes
!     with command MPI_BCAST
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBROADC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(201)
      CALL MPI_BCAST ( IPTR, ILEN, ITYPE, MASTER-1,
     &                 WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS = INTSTR(IERR)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(IF:IL)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF
!TIMG      CALL SWTSTO(201)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWGATHER ( IOPTR, IOLEN, IIPTR, IILEN, ITYPE )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Gathers different amounts of data from each processor
!     to the master
!
!  3. Method
!
!     Wrapper for MPI_GATHERV
!
!  4. Argument variables
!
!     IILEN       length of input array
!     IIPTR       pointer to first element of input array (local)
!     IOLEN       length of output array
!     IOPTR       pointer to first element of output array (global)
!     ITYPE       type of data
!
      INTEGER IILEN, IIPTR, IOLEN, IOPTR, ITYPE
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     ICOUNT:     array specifying array size of data received
!                 from each processor
!     IDSPLC:     array specifying the starting address of the
!                 incoming data from each processor, relative
!                 to the global array
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER                 I, IENT, IERR, IF, IL
      INTEGER, ALLOCATABLE :: ICOUNT(:), IDSPLC(:)
      CHARACTER*20            INTSTR, CHARS
      CHARACTER*80            MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_GATHER       Gathers data from all nodes to the master
!     MPI_GATHERV      Gathers different amounts of data from
!                      all nodes to the master
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     gather the array sizes to the master
!
!     check whether enough space has been allocated
!     for gathered data
!
!     calculate starting address of each local array
!     with respect to the global array
!
!     gather different amounts of data from each processor
!     to the master
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWGATHER')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF (IAMMASTER) THEN
         ALLOCATE(ICOUNT(0:NPROC-1))
         ALLOCATE(IDSPLC(0:NPROC-1))
      END IF

!     --- gather the array sizes to the master

      CALL MPI_GATHER( IILEN, 1, SWINT, ICOUNT, 1, SWINT,
     &                 MASTER-1, WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS = INTSTR(IERR)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(IF:IL)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

!     --- check whether enough space has been allocated
!         for gathered data

      IF (IAMMASTER) THEN
         IF ( SUM(ICOUNT).GT.IOLEN ) THEN
            CALL MSGERR(4,
     &                  'Not enough space allocated for gathered data')
            RETURN
         END IF
      END IF

!     --- calculate starting address of each local array
!         with respect to the global array

      IF (IAMMASTER) THEN
         IDSPLC(0) = 0
         DO I = 1, NPROC-1
            IDSPLC(I) = ICOUNT(I-1) + IDSPLC(I-1)
         END DO
      END IF

!     --- gather different amounts of data from each processor
!         to the master

      CALL MPI_GATHERV( IIPTR, IILEN, ITYPE, IOPTR, ICOUNT, IDSPLC,
     &                  ITYPE, MASTER-1, WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS = INTSTR(IERR)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(IF:IL)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF

      IF (IAMMASTER) DEALLOCATE(ICOUNT,IDSPLC)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCE ( IPTR, ILEN, ITYPE, ITYPRD )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.96, Dec. 08: call SWREDUCI/R instead of passing startaddress of the array
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on
!     array (IPTR) of type ITYPE to collect values from
!     all processes
!
!  4. Argument variables
!
!     ILEN        length of array to be collect
!     IPTR        pointer to first element of array to be collect
!     ITYPE       type of data
!     ITYPRD      type of reduction
!
      INTEGER IPTR, ILEN, ITYPE, ITYPRD
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!
      INTEGER IENT
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Performs a global reduction of data across all nodes
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCE')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN
!
!     --- actual reduction of field array based on its type
      IF ( ITYPE.EQ.SWINT ) THEN
         CALL SWREDUCI ( IPTR, ILEN, ITYPRD )
      ELSE IF ( ITYPE.EQ.SWREAL ) THEN
         CALL SWREDUCR ( IPTR, ILEN, ITYPRD )
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCI ( IARR, ILEN, ITYPRD )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4
      USE M_PARALL
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on integer
!     array IARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     IARR        integer array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      INTEGER IARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ITEMP :     temporary array to store collected data
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      INTEGER, ALLOCATABLE :: ITEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_ALLREDUCE    Combines values from all processes and
!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Performs a global reduction of integers across all nodes
!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCI')

      ALLOCATE(ITEMP(ILEN))

!TIMG      CALL SWTSTA(202)
      CALL MPI_ALLREDUCE ( IARR, ITEMP, ILEN, SWINT,
     &                     ITYPRD, WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS = INTSTR(IERR)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(IF:IL)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF
      IARR = ITEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(ITEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCR ( ARR, ILEN, ITYPRD )
!
!****************************************************************
!
      !USE MPI
      USE OCPCOMM4
      USE M_PARALL
!
      IMPLICIT NONE
      include 'mpif.h'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on real
!     array ARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     ARR         real array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      REAL    ARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!     TEMP  :     temporary array to store collected data
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      REAL, ALLOCATABLE :: TEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MPI_ALLREDUCE    Combines values from all processes and
!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Performs a global reduction of reals across all nodes
!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCR')

      ALLOCATE(TEMP(ILEN))

!TIMG      CALL SWTSTA(202)
      CALL MPI_ALLREDUCE ( ARR, TEMP, ILEN, SWREAL,
     &                     ITYPRD, WAV_COMM_WORLD, IERR )
      IF ( IERR.NE.MPI_SUCCESS ) THEN
         CHARS = INTSTR(IERR)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'MPI produces some internal error - '//
     &            'return code is '//CHARS(IF:IL)
         CALL MSGERR ( 4, MSGSTR )
         RETURN
      END IF
      ARR = TEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(TEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSTRIP ( IPOWN, IDIR, NPART, IWORK, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Performs a stripwise partitioning with straight interfaces
!
!  3. Method
!
!     Each active point in a row/column will be assign to a part
!     according to its number and size (stored in IWORK).
!     The remaining points in the row/column will be assign
!     to the same part.
!
!  4. Argument variables
!
!     IDIR        direction of cutting
!                 1 = row
!                 2 = column
!     IPOWN       array giving the subdomain number of each gridpoint
!     IWORK       work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!     NPART       number of parts to be created
!
      INTEGER   IDIR, MXC, MYC, NPART
      INTEGER   IPOWN(*)
      INTEGER*8 IWORK(2,*)
!
!  6. Local variables
!
!     IC    :     index of (IX,IY)-point
!     ICC   :     index of (IX,IY)-point
!     IENT  :     number of entries
!     INCX  :     increment for adressing: 1 for x-dir, MXC for y-dir
!     INCY  :     increment for adressing: MXC for x-dir, 1 for y-dir
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IYY   :     index in y-direction
!     IPART :     a part counter
!     MXCI  :     maximum counter of gridpoints in x/y-direction
!     MYCI  :     maximum counter of gridpoints in y/x-direction
!     NCURPT:     number of currently assigned points to a created part
!     NPREM :     number of remaining points in a row/column
!
      INTEGER IC, ICC, IENT, INCX, INCY, IX, IY, IYY, IPART,
     &        MXCI, MYCI, NCURPT, NPREM
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWPARTIT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     depending on cutting direction, determine indirect addressing
!     create first empty part
!     for all active points do
!         assign this point to the created part
!         if size of created part has been reached
!            determine remaining active points in the current column
!            if no remaining points, create next empty part
!            else remaining points belong to the current part
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSTRIP')

!     --- depending on cutting direction, determine indirect addressing
!         for array IPOWN

      IF ( IDIR.EQ.1 ) THEN
         MXCI = MYC
         MYCI = MXC
         INCX = MXC
         INCY = 1
      ELSE IF ( IDIR.EQ.2 ) THEN
         MXCI = MXC
         MYCI = MYC
         INCX = 1
         INCY = MXC
      END IF

!     --- create first empty part

      IPART  = 1
      NCURPT = 0

!     --- for all active points do

      DO IX = 1, MXCI
         DO IY = 1, MYCI

            IC = IX*INCX + IY*INCY - MXC

            IF ( IPOWN(IC).EQ.1 ) THEN

!              --- assign this point to the created part

               IPOWN(IC) = IWORK(1,IPART)
               NCURPT    = NCURPT + 1

!              --- if size of created part has been reached

               IF ( NCURPT.GE.IWORK(2,IPART) ) THEN

!                 --- determine remaining active points in the
!                     current column

                  NPREM = 0
                  DO IYY = IY+1, MYCI
                     ICC = IX*INCX + IYY*INCY - MXC
                     IF (IPOWN(ICC).EQ.1) NPREM = NPREM +1
                  END DO

                  IF ( NPREM.EQ.0 ) THEN

!                    --- if no remaining points, create next empty part

                     IPART  = IPART + 1
                     NCURPT = 0

                  ELSE

!                    --- else remaining points belong to the current part

                     IWORK(2,IPART  ) = IWORK(2,IPART  ) + NPREM
                     IWORK(2,IPART+1) = IWORK(2,IPART+1) - NPREM

                  END IF

               END IF

            END IF

         END DO
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWPARTIT ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Sep. 04: determines load per processor based on speed
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Carries out the partitioning of the SWAN computational grid
!
!  3. Method
!
!     Based on stripwise partitioning
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  6. Local variables
!
!     I     :     loop counter
!     ICNT  :     auxiliary integer to count weights
!     IDIR  :     direction of cutting
!                 1 = row
!                 2 = column
!     IENT  :     number of entries
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IWORK :     work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     NACTP :     total number of active gridpoints
!     NPCUM :     cumulative number of gridpoints
!
      INTEGER   I, IDIR, IENT, IX, IY
      INTEGER*8 ICNT, NACTP, NPCUM
      INTEGER*8 IWORK(2,NPROC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWSTRIP          Performs a stripwise partitioning with straight
!                      interfaces
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWDECOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!     determine direction of cutting
!     determine number of active points
!     determine numbers and sizes of parts to be created
!     partition grid
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWPARTIT')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- determine direction of cutting

      IF ( MXC.GT.MYC ) THEN
         IDIR = 2
      ELSE
         IDIR = 1
      END IF

!     --- determine number of active points and
!         set IPOWN to 1 in these points

      NACTP = 0
      DO IX = 1, MXC
         DO IY = 1, MYC
            IF ( KGRPGL(IX,IY).NE.1 ) THEN
               IPOWN(IX,IY) = 1
               NACTP        = NACTP + 1
            END IF
         END DO
      END DO

!     --- determine numbers and sizes of parts to be created

      NPCUM = 0
      ICNT  = 0
      DO I = 1, NPROC
         ICNT       = ICNT + IWEIG(I)
         IWORK(1,I) = I
         IWORK(2,I) = (NACTP*ICNT)/SUM(IWEIG) - NPCUM
         NPCUM      = (NACTP*ICNT)/SUM(IWEIG)
      END DO
      DEALLOCATE(IWEIG)

!     --- partition grid

      CALL SWSTRIP ( IPOWN, IDIR, NPROC, IWORK, MXC, MYC )

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBLADM ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Jul. 04: determine global bounds in subdomains
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     For the present node, carries out the block administration
!     and determines array bounds with respect to global grid
!
!  3. Method
!
!     Based on domain decomposition, the interface sizes are
!     determined that is needed for the setup of block
!     administration stored as IBLKAD
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     IC    :     index of (IX,IY)-point
!     ICOFF :     offset of IC-index
!     ICRECV:     array containing positions of unknowns
!                 to be received from neighbour
!     ICSEND:     array containing positions of unknowns
!                 to be sent to neighbour
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     INB   :     neighbour counter
!     ISTART:     startaddress for each size interface in array IBLKAD
!     IX    :     index in x-direction
!     IXOFF :     offset in x-direction
!     IY    :     index in y-direction
!     IYOFF :     offset in y-direction
!     IWORK :     array used to determine interface sizes
!                   IWORK(1,i) = number of the i-th neighbour
!                   IWORK(2,i) = position of the i-th neighbour with
!                                respect to present subdomain
!                                (resp. top, bottom, right, left)
!                   IWORK(3,i) = size of interface to i-th neighbour
!     JOFFS :     offsets at which a point of a neigbhour domain can be found
!     MSGSTR:     string to pass message to call MSGERR
!     MXSIZ :     size of present subdomain in x-direction
!     MYSIZ :     size of present subdomain in y-direction
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!
      INTEGER      I, IC, ICOFF, IDOM, IENT, IF, IL, INB, ISTART,
     &             IX, IXOFF, IY, IYOFF, JOFFS(2,4), MXSIZ, MYSIZ,
     &             NNEIGH, NOVLU
      INTEGER      IWORK(3,NPROC),
     &             ICRECV(NPROC,MAX(MXC,MYC)),
     &             ICSEND(NPROC,MAX(MXC,MYC))
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWDECOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     intialize offsets to be used in searching for interfaces
!     determine enclosing box of present subdomain
!     if subdomain appears to be empty
!        give warning and set empty bounding box
!     else
!        extend enclosing box to include halo area
!     localize global bounds in present subdomain
!     determine size of enclosing box
!     determine interface sizes:
!
!        loop over global grid
!           if point belongs to this part
!              for each of the four sizes
!                  if a neighbouring subdomain is found there
!                     find it in the list of neighbours
!                     if not yet in the list, add it
!                     store position of neighbour
!                     update number of overlapping unknowns
!
!     store block administration
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBLADM')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- intialize offsets to be used in searching for interfaces

      JOFFS = RESHAPE((/0,1,0,-1,1,0,-1,0/), (/2,4/))

!     --- determine enclosing box of present subdomain

      IF ( MXC.GT.MYC ) THEN
         MXF = MXC+1
         MXL = 0
      ELSE
         MYF = MYC+1
         MYL = 0
      END IF

      DO IX = 1, MXC
         DO IY = 1, MYC

            IF( IPOWN(IX,IY).EQ.INODE ) THEN

               MXF = MIN(IX,MXF)
               MYF = MIN(IY,MYF)
               MXL = MAX(IX,MXL)
               MYL = MAX(IY,MYL)

            END IF

         END DO
      END DO

!     --- if subdomain appears to be empty

      IF ( MXF.GT.MXL .OR. MYF.GT.MYL ) THEN

!        --- give warning and set empty bounding box

         CHARS = INTSTR(INODE)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'Empty subdomain is detected - '//
     &            ' node number is '//CHARS(IF:IL)
         CALL MSGERR ( 2, MSGSTR )

         MXF = 1
         MYF = 1
         MXL = 0
         MYL = 0

      ELSE

!        --- extend enclosing box to include halo area

         MXF = MAX(1  ,MXF-IHALOX)
         MYF = MAX(1  ,MYF-IHALOY)
         MXL = MIN(MXC,MXL+IHALOX)
         MYL = MIN(MYC,MYL+IHALOY)

      END IF

!     --- localize global bounds in present subdomain

      IF ( MXCGL.GT.MYCGL ) THEN
         LMXF = MXF.EQ.1     .AND. INODE.EQ.1
         LMXL = MXL.EQ.MXCGL .AND. INODE.EQ.NPROC
         LMYF = MYF.EQ.1
         LMYL = MYL.EQ.MYCGL
      ELSE
         LMXF = MXF.EQ.1
         LMXL = MXL.EQ.MXCGL
         LMYF = MYF.EQ.1     .AND. INODE.EQ.1
         LMYL = MYL.EQ.MYCGL .AND. INODE.EQ.NPROC
      END IF

!     --- determine size of enclosing box

      MXSIZ = MXL - MXF + 1
      MYSIZ = MYL - MYF + 1

      IWORK  = 0
      ICRECV = 0
      ICSEND = 0

!     --- determine interface sizes

      DO IX = 1, MXC
         DO IY = 1, MYC

!           --- if point belongs to this part

            IF ( IPOWN(IX,IY).EQ.INODE ) THEN

!              --- for each of the four sizes

               DO I = 1, 4

                  IXOFF = JOFFS(1,I)
                  IYOFF = JOFFS(2,I)

!                 --- if a neighbouring subdomain is found there

                  IF ( (IX+IXOFF).GT.0.AND.(IX+IXOFF).LE.MXC.AND.
     &                 (IY+IYOFF).GT.0.AND.(IY+IYOFF).LE.MYC ) THEN

                     IF ( IPOWN(IX+IXOFF,IY+IYOFF).NE.0.AND.
     &                    IPOWN(IX+IXOFF,IY+IYOFF).NE.INODE ) THEN

                        IC    = (IY      -MYF)*MXSIZ + (IX      -MXF+1)
                        ICOFF = (IY+IYOFF-MYF)*MXSIZ + (IX+IXOFF-MXF+1)

!                       --- find it in the list of neighbours

                        IDOM = IPOWN(IX+IXOFF,IY+IYOFF)

                        INB = 1
  100                   IF ( INB.LE.NPROC .AND.
     &                       IWORK(1,INB).NE.IDOM .AND.
     &                       IWORK(1,INB).NE.0 )  THEN
                           INB = INB + 1
                           GOTO 100
                        END IF

                        IF ( INB.GT.NPROC ) THEN
                          CALL MSGERR (4,'Found more neighbours than '//
     &                                 'subdomains in the partitioning')
                          RETURN
                        END IF

!                       --- if not yet in the list, add it

                        IF ( IWORK(1,INB).EQ.0 ) IWORK(1,INB) = IDOM

!                       --- store position of neighbour with respect to
!                           present subdomain

                        IWORK(2,INB) = I

!                       --- update number of overlapping unknowns

                        IWORK(3,INB) = IWORK(3,INB) + 1

                        ICSEND(INB,IWORK(3,INB)) = IC
                        ICRECV(INB,IWORK(3,INB)) = ICOFF

                     END IF

                  END IF

               END DO

            END IF

         END DO
      END DO

!     --- store block administration

      NNEIGH    = COUNT(IWORK(1,:)>0)
      IBLKAD(1) = NNEIGH
      ISTART    = 3*NNEIGH+2
      DO INB = 1, NNEIGH
         IBLKAD(3*INB-1) = IWORK(1,INB)
         IBLKAD(3*INB  ) = IWORK(2,INB)
         IBLKAD(3*INB+1) = ISTART
         NOVLU           = IWORK(3,INB)
         IBLKAD(ISTART)  = NOVLU
         DO I = 1, NOVLU
            IBLKAD(ISTART      +I) = ICSEND(INB,I)
            IBLKAD(ISTART+NOVLU+I) = ICRECV(INB,I)
         END DO
         ISTART = ISTART + 2*NOVLU+1
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWDECOMP
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Carries out domain decomposition meant for
!     distributed-memory approach
!
!  3. Method
!
!     First, carry out the partitioning of the
!     SWAN computational grid and then do the
!     block administration
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!     IX    :     loop counter
!     IY    :     loop counter
!     IPOWN :     array giving the subdomain number of each gridpoint
!
      INTEGER IENT, IX, IY
      INTEGER, ALLOCATABLE :: IPOWN(:,:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWBLADM          Carries out the block administration
!     SWPARTIT         Carries out the partitioning of the SWAN
!                      computational grid
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     allocate and initialize array for block administration
!     store the original values of MXC, MYC and MCGRD
!     if not parallel, return
!     carry out the partitioning of computational grid
!     carry out the block administration
!     compute MXC, MYC and MCGRD for each subdomain
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWDECOMP')

!     --- allocate and initialize array for block administration

      IF (.NOT.ALLOCATED(IBLKAD)) ALLOCATE(IBLKAD(41+20*MAX(MXC,MYC)))
      IBLKAD = 0

!     --- store the original values of MXC, MYC and MCGRD of
!         global computational grid

      MXCGL   = MXC
      MYCGL   = MYC
      MCGRDGL = MCGRD
      MXF     = 1
      MXL     = MXC
      MYF     = 1
      MYL     = MYC
      LMXF    = MXF.EQ.1                                                  40.41
      LMXL    = MXL.EQ.MXCGL                                              40.41
      LMYF    = MYF.EQ.1                                                  40.41
      LMYL    = MYL.EQ.MYCGL                                              40.41

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- carry out the partitioning of the SWAN
!         computational grid

      ALLOCATE(IPOWN(MXC,MYC))
      IPOWN = 0
      CALL SWPARTIT( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- carry out the block administration

      CALL SWBLADM( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- compute MXC, MYC and MCGRD for each subdomain

      MXC = MXL - MXF + 1
      MYC = MYL - MYF + 1

      MCGRD = 1
      DO IX = MXF, MXL
         DO IY = MYF, MYL
            IF ( KGRPGL(IX,IY).NE.1 ) MCGRD = MCGRD + 1
         END DO
      END DO

      DEALLOCATE(IPOWN)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXCHG ( FIELD, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Updates geographical field array through exchanging
!     values between neighbouring subdomains
!
!  3. Method
!
!     Made use of MPI by means of SWSENDNB and SWRECVNB
!     and also block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     FIELD       geographical field array for which 'halo' values must
!                 be copied from neighbouring subdomains
!     KGRPNT      indirect addressing for grid points
!
      INTEGER KGRPNT(MXC*MYC)
      REAL    FIELD(MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     ISTART:     pointer in array IBLKAD
!     ITAG  :     message tag for sending and receiving
!     K     :     loop counter
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     WORK  :     work array to store data to be sent to or
!                 received from neighbour
!
      INTEGER IDOM, IENT, INB, ISTART, ITAG, K, NNEIGH, NOVLU
      REAL    WORK(MAX(MXC,MYC))
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!     SWSENDNB         Data is sent to a neighbour
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get subdomain number, pointer and size
!        store data to be sent in array WORK
!        send array WORK
!
!     for all neighbouring subdomains do
!        get subdomain number, pointer and size
!        receive next array and store in WORK
!        store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWEXCHG')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(203)

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get subdomain number, pointer and size

         IDOM   = IBLKAD(3*INB-1)
         ISTART = IBLKAD(3*INB+1)
         NOVLU  = IBLKAD(ISTART)

!        --- store data to be sent in array WORK

         DO K = 1, NOVLU
            WORK(K) = FIELD(KGRPNT(IBLKAD(ISTART+K)))
         END DO

!        --- send array WORK

         ITAG = 2
         CALL SWSENDNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
         IF (STPNOW()) RETURN

      END DO

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get subdomain number, pointer and size

         IDOM   = IBLKAD(3*INB-1)
         ISTART = IBLKAD(3*INB+1)
         NOVLU  = IBLKAD(ISTART)

!        --- receive next array and store in WORK

         ITAG  = 2
         CALL SWRECVNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
         IF (STPNOW()) RETURN

!        --- store the received data

         DO K = 1, NOVLU
            FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+K))) = WORK(K)
         END DO

      END DO

!TIMG      CALL SWTSTO(203)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVAC ( AC2, IS, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Receives action density from neighbouring subdomains
!     depending on sweep direction
!
!  3. Method
!
!     Use of SWRECVNB and block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     AC2         action density
!     IS          start index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     SWPDIR      sweep direction
!
      INTEGER IS, J, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    AC2(MDC,MSC,MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     IPR   :     array containing positions of neighbours from
!                 which data is to be received
!     ITAG  :     message tag for sending and receiving
!     NNEIGH:     number of neighbouring subdomains
!     WORK  :     work array to store data to received from neighbour
!
      INTEGER IDOM, IENT, INB, IPNB, ITAG, NNEIGH
      INTEGER IPR(2,4)
      REAL    WORK(MDC,MSC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if position corresponds to sweep selection
!           get subdomain number
!           receive next array and store in WORK
!           store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVAC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IPR = RESHAPE((/2,4,2,3,1,3,1,4/), (/2,4/))

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if position corresponds to sweep selection

         IF ( IPNB.EQ.IPR(1,SWPDIR) .OR. IPNB.EQ.IPR(2,SWPDIR) ) THEN

!           --- get subdomain number

            IDOM   = IBLKAD(3*INB-1)

!           --- receive next array and store in WORK

            ITAG  = 2
            CALL SWRECVNB ( WORK, MDC*MSC, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

!           --- store the received data

            AC2(:,:,KGRPNT(IS,J)) = WORK(:,:)

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDAC ( AC2, IE, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Sends action density to neighbouring subdomains
!     depending on sweep direction
!
!  3. Method
!
!     Use of SWSENDNB and block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     AC2         action density
!     IE          end index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     SWPDIR      sweep direction
!
      INTEGER IE, J, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    AC2(MDC,MSC,MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     IPS   :     array containing positions of neighbours to
!                 which data is to be sent
!     ITAG  :     message tag for sending and receiving
!     NNEIGH:     number of neighbouring subdomains
!     WORK  :     work array to store data to be sent to neighbour
!
      INTEGER IDOM, IENT, INB, IPNB, ITAG, NNEIGH
      INTEGER IPS(2,4)
      REAL    WORK(MDC,MSC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWSENDNB         Data is sent to a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if position corresponds to sweep selection
!           get subdomain number
!           store data to be sent in array WORK
!           send array WORK
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDAC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IPS = RESHAPE((/1,3,1,4,2,4,2,3/), (/2,4/))

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if position corresponds to sweep selection

         IF ( IPNB.EQ.IPS(1,SWPDIR) .OR. IPNB.EQ.IPS(2,SWPDIR) ) THEN

!           --- get subdomain number

            IDOM   = IBLKAD(3*INB-1)

!           --- store data to be sent in array WORK

            WORK(:,:) = AC2(:,:,KGRPNT(IE,J))

!           --- send array WORK

            ITAG = 2
            CALL SWSENDNB ( WORK, MDC*MSC, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLLECT ( FIELDGL, FIELD, FULL )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_GENARR                                                        40.31
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Mar. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: extended to full arrays
!
!  2. Purpose
!
!     Collects geographical FIELD arrays from all nodes
!
!  3. Method
!
!     Made use of MPI by means of SWGATHER
!
!  4. Argument variables
!
!     FIELD       geographical field array in own subdomain
!     FIELDGL     global geographical field array gathered from all nodes
!     FULL        if true, full arrays are handled otherwise 1-D compact
!                 arrays are handled
!
      REAL    FIELD(*), FIELDGL(*)
      LOGICAL FULL
!
!  6. Local variables
!
!     FLDC  :     auxiliary array for collecting data
!     IARRC :     auxiliary array for collecting grid indices and counter
!     IARRL :     auxiliary array containing grid indices and counter
!     IENT  :     number of entries
!     ILEN  :     integer indicating length of an array
!     ILEN2 :     integer indicating length of another array
!     INDX  :     pointer in array
!     INDXC :     pointer in collected array
!     IOFF1 :     offset
!     IOFF2 :     another offset
!     IP    :     node number
!     IX    :     loop counter
!     IY    :     loop counter
!     KGRPTC:     auxiliary array for collecting indirect
!                 addressing for grid points
!     MXFGL :     first index w.r.t. global grid in x-direction
!     MXLGL :     last index w.r.t. global grid in x-direction
!     MYFGL :     first index w.r.t. global grid in y-direction
!     MYLGL :     last index w.r.t. global grid in y-direction
!
      INTEGER IENT, ILEN, ILEN2, INDX, INDXC, IOFF1, IOFF2, IP, IX, IY,
     &        MXFGL, MXLGL, MYFGL, MYLGL
      INTEGER IARRL(5), IARRC(5,0:NPROC-1)

      INTEGER, ALLOCATABLE :: KGRPTC(:)
      REAL,    ALLOCATABLE :: FLDC(:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWGATHER         Gathers different amounts of data from all nodes
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if sequential run
!        just make a copy
!     else
!        gather necessary arrays
!        copy gathered data to global array in appropriate manner
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLLECT')

      IF (.NOT.PARLL) THEN

!        --- in case of sequential run, just make a copy

         IF (FULL) THEN
            DO INDX = 1, MXC*MYC
               FIELDGL(INDX) = FIELD(INDX)
            END DO
         ELSE
            DO INDX = 1, MCGRD
               FIELDGL(INDX) = FIELD(INDX)
            END DO
         END IF

      ELSE

!        --- gather necessary arrays

         IARRL(1) = MXF
         IARRL(2) = MXL
         IARRL(3) = MYF
         IARRL(4) = MYL
         IARRL(5) = MCGRD
         CALL SWGATHER (IARRC, 5*NPROC, IARRL, 5, SWINT )
         IF (STPNOW()) RETURN

         IF (.NOT.FULL) THEN
            IF (IAMMASTER) THEN
               ILEN = MXCGL*MYCGL +
     &                       4*MAX(IHALOX,IHALOY)*NPROC*MAX(MXCGL,MYCGL)
               ALLOCATE(KGRPTC(ILEN))
            END IF
            CALL SWGATHER ( KGRPTC, ILEN, KGRPNT, MXC*MYC, SWINT )
            IF (STPNOW()) RETURN
         ELSE
            IF (IAMMASTER) ALLOCATE(KGRPTC(0))
         END IF

         IF (IAMMASTER) THEN
            IF (FULL) THEN
               ILEN = MXCGL*MYCGL +
     &                       4*MAX(IHALOX,IHALOY)*NPROC*MAX(MXCGL,MYCGL)
            ELSE
               ILEN = MCGRDGL +
     &                       4*MAX(IHALOX,IHALOY)*NPROC*MAX(MXCGL,MYCGL)
            END IF
            ALLOCATE(FLDC(ILEN))
         END IF
         IF (FULL) THEN
            ILEN2 = MXC*MYC
         ELSE
            ILEN2 = MCGRD
         END IF
         CALL SWGATHER ( FLDC, ILEN, FIELD, ILEN2, SWREAL )
         IF (STPNOW()) RETURN

!        --- copy gathered data to global array in appropriate manner

         IF (IAMMASTER) THEN

            IOFF1 = 0
            IOFF2 = 0

            DO IP = 0, NPROC-1
               IF ( MXCGL.GT.MYCGL ) THEN
                  IF ( IARRC(1,IP).EQ.1 .AND. IP.EQ.0 ) THEN
                     MXFGL = 1
                  ELSE
                     MXFGL = IARRC(1,IP) + IHALOX
                  END IF
                  IF ( IARRC(3,IP).EQ.1 ) THEN
                     MYFGL = 1
                  ELSE
                     MYFGL = IARRC(3,IP) + IHALOY
                  END IF
                  IF ( IARRC(2,IP).EQ.MXCGL .AND. IP.EQ.NPROC-1 ) THEN
                     MXLGL = MXCGL
                  ELSE
                     MXLGL = IARRC(2,IP) - IHALOX
                  END IF
                  IF ( IARRC(4,IP).EQ.MYCGL ) THEN
                     MYLGL = MYCGL
                  ELSE
                     MYLGL = IARRC(4,IP) - IHALOY
                  END IF
               ELSE
                  IF ( IARRC(1,IP).EQ.1 ) THEN
                     MXFGL = 1
                  ELSE
                     MXFGL = IARRC(1,IP) + IHALOX
                  END IF
                  IF ( IARRC(3,IP).EQ.1 .AND. IP.EQ.0 ) THEN
                     MYFGL = 1
                  ELSE
                     MYFGL = IARRC(3,IP) + IHALOY
                  END IF
                  IF ( IARRC(2,IP).EQ.MXCGL ) THEN
                     MXLGL = MXCGL
                  ELSE
                     MXLGL = IARRC(2,IP) - IHALOX
                  END IF
                  IF ( IARRC(4,IP).EQ.MYCGL .AND. IP.EQ.NPROC-1 ) THEN
                     MYLGL = MYCGL
                  ELSE
                     MYLGL = IARRC(4,IP) - IHALOY
                  END IF
               END IF

               ILEN = IARRC(2,IP)-IARRC(1,IP)+1

               IF (FULL) THEN
                  DO IX = MXFGL, MXLGL
                     DO IY = MYFGL, MYLGL
                        INDX  = (IY-1)*MXCGL+IX
                        INDXC = (IY-IARRC(3,IP))*ILEN+IX
     &                             -IARRC(1,IP)+1+IOFF2
                        FIELDGL(INDX)= FLDC(INDXC)
                     END DO
                  END DO
               ELSE
                  DO IX = MXFGL, MXLGL
                     DO IY = MYFGL, MYLGL
                        INDX  = KGRPGL(IX,IY)
                        INDXC = KGRPTC((IY-IARRC(3,IP))*ILEN+IX
     &                                    -IARRC(1,IP)+1+IOFF2)
                        FIELDGL(INDX) = FLDC(INDXC+IOFF1)
                     END DO
                  END DO
               END IF

               IOFF1 = IOFF1 + IARRC(5,IP)
               IOFF2 = IOFF2 + (IARRC(2,IP)-IARRC(1,IP)+1)*
     &                         (IARRC(4,IP)-IARRC(3,IP)+1)

            END DO

         END IF

         IF (IAMMASTER) DEALLOCATE(KGRPTC,FLDC)

      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLOUT ( OURQT, BLKND )                                40.51
!
!****************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.51
      USE SWCOMM3                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: re-design output process in parallel mode
!
!  2. Purpose
!
!     Collects output results
!
!  3. Method
!
!     Read individual process output files containing tables,
!     spectral outputs and block data and write them to
!     generic output files in appropriate manner
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain        40.41
!     OURQT       array indicating at what time requested output          40.51
!                 is processed                                            40.51
!
      REAL    BLKND(MXCGL,MYCGL), OURQT(MAX_OUTP_REQ)                     40.51
!
!  6. Local variables
!
!     CORQ  :     current item in list of request outputs
!     CROSS :     auxiliary logical array
!     CUOPS :     current item in list of point sets
!     DIF   :     difference between end and actual times
!     DTTIWR:     to write time string
!     IC    :     loop variable
!     IENT  :     number of entries
!     IP    :     loop variable
!     IRQ   :     request number
!     IT    :     time step counter
!     IT0   :     integer indicating first step of simulation
!     IT1   :     integer indicating last step of simulation
!     ITMP1 :     auxiliary integer
!     ITMP2 :     auxiliary integer
!     ITMP3 :     auxiliary integer
!     ITMP4 :     auxiliary integer
!     ITMP5 :     auxiliary integer
!     ITMP6 :     auxiliary integer
!     IUNIT :     counter for file unit numbers
!     MIP   :     total number of output points
!     MXK   :     number of points in x-direction of output frame
!     MYK   :     number of points in y-direction of output frame
!     OPENED:     logical whether a file is open or not
!     PSTYPE:     type of point set
!     RTMP1 :     auxiliary real
!     RTMP2 :     auxiliary real
!     RTMP3 :     auxiliary real
!     RTMP4 :     auxiliary real
!     RTYPE :     type of request
!     SNAMPF:     name of plot frame
!     TNEXT :     time of next requested output
!     XC    :     computational grid x-coordinate of output point
!     XP    :     user x-coordinate of output point
!     YC    :     computational grid y-coordinate of output point
!     YP    :     user y-coordinate of output point
!
      INTEGER   IC, IENT, IP, IRQ, IT, IT0, IT1, IUNIT, MIP, MXK, MYK
      INTEGER   ITMP1, ITMP2, ITMP3, ITMP4, ITMP5, ITMP6
      REAL      DIF, TNEXT
      REAL      RTMP1, RTMP2, RTMP3, RTMP4
      REAL, ALLOCATABLE :: XC(:), YC(:), XP(:), YP(:)                     40.51
      LOGICAL   OPENED
      LOGICAL, ALLOCATABLE :: CROSS(:,:)                                  40.86
      CHARACTER PSTYPE*1, RTYPE*4, SNAMPF*8, DTTIWR*18
      TYPE(ORQDAT), POINTER :: CORQ
      TYPE(OPSDAT), POINTER :: CUOPS
!
!  8. Subroutines used
!
!     EQREAL           Logical comparing two reals
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWCOLBLK         Collects block output
!     SWCOLSPC         Collects spectral output
!     SWCOLTAB         Collects table ouput
!     SWOEXC           Computes coordinates of output points              40.51
!
      LOGICAL   EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     do for all COMPUTE commands
!        do for all time steps
!           do for all output requests
!              processing of output instructions necessary for collection
!              check time of output action
!              compute coordinates of output points
!              correct problem coordinates with offset values
!              rewrite table output by means of collection of output
!              rewrite spectral output by means of collection of output
!              rewrite block output by means of collection of output
!     close all files and delete process files
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLOUT')

      IF ( NREOQ.EQ.0 ) RETURN

!     --- do for all COMPUTE commands

      DO IC = 1, NCOMPT

         NSTATC = NINT(RCOMPT(IC,1))
         IF ( NSTATC.EQ.1 ) THEN
            IT0 = 0
         ELSE
            IT0 = 1
         END IF
         IT1   = NINT(RCOMPT(IC,2))
         TFINC = RCOMPT(IC,3)
         TINIC = RCOMPT(IC,4)
         DT    = RCOMPT(IC,5)
         TIMCO = TINIC

!        --- do for all time steps

         DO IT = IT0, IT1

            IF (NSTATM.GT.0) CHTIME = DTTIWR(ITMOPT, TIMCO)

!           --- do for all output requests

            CORQ => FORQ
            DO 100 IRQ = 1, NREOQ

!              --- processing of output instructions necessary
!                  for collection

!              --- check time of output action

               DIF = TFINC - TIMCO
               IF ( IT.EQ.IT0 .AND. IC.EQ.1 ) THEN
                  CORQ%OQR(1) = OURQT(IRQ)                                40.51
               END IF
               TNEXT = CORQ%OQR(1)
               IF ( ABS(DIF).LT.0.5*DT .AND. CORQ%OQR(2).LT.0. ) THEN
                  CORQ%OQR(1) = TIMCO
               ELSE IF ( CORQ%OQR(2).GT.0. .AND. TIMCO.GE.TNEXT ) THEN
                  CORQ%OQR(1) = TNEXT + CORQ%OQR(2)
               ELSE
                  GOTO 50
               END IF

               RTYPE  = CORQ%RQTYPE
               SNAMPF = CORQ%PSNAME

               CUOPS => FOPS
               DO
                 IF (CUOPS%PSNAME.EQ.SNAMPF) EXIT
                 IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) GOTO 50
                 CUOPS => CUOPS%NEXTOPS
               END DO
               PSTYPE = CUOPS%PSTYPE

               IF ( PSTYPE.EQ.'F' .OR. PSTYPE.EQ.'H' ) THEN
                  MXK = CUOPS%OPI(1)
                  MYK = CUOPS%OPI(2)
                  MIP = MXK * MYK
               ELSE IF ( PSTYPE.EQ.'C' .OR. PSTYPE.EQ.'P' .OR.
     &                   PSTYPE.EQ.'N' ) THEN
                  MXK = 0
                  MYK = 0
                  MIP = CUOPS%MIP
               ELSE IF ( PSTYPE.EQ.'U' ) THEN
                  MIP = CUOPS%MIP
                  MXK = MIP
                  MYK = 1
               END IF

               IF (.NOT.ALLOCATED(XC)) ALLOCATE(XC(MIP))
               IF (.NOT.ALLOCATED(YC)) ALLOCATE(YC(MIP))
               IF (.NOT.ALLOCATED(XP)) ALLOCATE(XP(MIP))
               IF (.NOT.ALLOCATED(YP)) ALLOCATE(YP(MIP))
               IF (.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(4,MIP))          40.86

!              --- compute coordinates of output points                   40.51

               ITMP1  = MXC
               ITMP2  = MYC
               ITMP3  = MCGRD
               ITMP4  = NGRBND
               ITMP5  = MXF
               ITMP6  = MYF
               RTMP1  = XCLMIN
               RTMP2  = XCLMAX
               RTMP3  = YCLMIN
               RTMP4  = YCLMAX
               MXC    = MXCGL
               MYC    = MYCGL
               MCGRD  = MCGRDGL
               NGRBND = NGRBGL
               MXF    = 1
               MYF    = 1
               XCLMIN = XCGMIN
               XCLMAX = XCGMAX
               YCLMIN = YCGMIN
               YCLMAX = YCGMAX
               CALL SWOEXC (PSTYPE              ,
     &                      CUOPS%OPI           ,CUOPS%OPR           ,
     &                      CUOPS%XP            ,CUOPS%YP            ,
     &                      MIP                 ,XP                  ,
     &                      YP                  ,XC                  ,
     &                      YC                  ,KGRPGL              ,
     &                      XGRDGL              ,YGRDGL              ,
     &                      CROSS               )                         40.86
               MXC    = ITMP1
               MYC    = ITMP2
               MCGRD  = ITMP3
               NGRBND = ITMP4
               MXF    = ITMP5
               MYF    = ITMP6
               XCLMIN = RTMP1
               XCLMAX = RTMP2
               YCLMIN = RTMP3
               YCLMAX = RTMP4

!              --- correct problem coordinates with offset values

               DO IP = 1, MIP
                  RTMP1 = XP(IP)
                  RTMP2 = YP(IP)
                  IF (.NOT.EQREAL(RTMP1,OVEXCV(1))) XP(IP)=RTMP1+XOFFS
                  IF (.NOT.EQREAL(RTMP2,OVEXCV(2))) YP(IP)=RTMP2+YOFFS
               END DO

!              --- rewrite table output by means of collection of
!                  output locations

               IF ( RTYPE(1:3).EQ.'TAB' ) THEN
                  CALL SWCOLTAB ( RTYPE, CORQ%OQI, CORQ%IVTYP, MIP, IRQ,
     &                            BLKND, XC, YC, XP, YP )                 40.51 40.41
                  IF (STPNOW()) RETURN
               END IF

!              --- rewrite spectral output by means of collection of
!                  output locations

               IF ( RTYPE(1:2).EQ.'SP' ) THEN
                  CALL SWCOLSPC ( RTYPE, CORQ%OQI, MIP, IRQ, BLKND,       40.51 40.41
     &                            XC, YC )                                40.51
                  IF (STPNOW()) RETURN
               END IF

!              --- rewrite block output by means of collection of process
!                  output data

               IF ( RTYPE(1:3).EQ.'BLK' ) THEN
                  CALL SWCOLBLK ( RTYPE, CORQ%OQI, CORQ%IVTYP, CORQ%FAC,
     &                            SNAMPF, MXK, MYK, IRQ, BLKND, XC, YC )  40.51 40.41
                  IF (STPNOW()) RETURN
               END IF

               IF (ALLOCATED(XP)) DEALLOCATE(XP)                          40.51
               IF (ALLOCATED(YP)) DEALLOCATE(YP)                          40.51
               IF (ALLOCATED(XC)) DEALLOCATE(XC)                          40.51
               IF (ALLOCATED(YC)) DEALLOCATE(YC)                          40.51

               IF (ALLOCATED(CROSS)) DEALLOCATE(CROSS)                    40.86

  50           CONTINUE
               CORQ => CORQ%NEXTORQ

 100        CONTINUE

            IF ( NSTATC.EQ.1.AND.IT.LT.IT1 ) TIMCO = TIMCO + DT

         END DO

      END DO

!     --- close all files and delete process files

      DO IUNIT = HIOPEN+1, HIOPEN+NREOQ
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE(IUNIT)
      END DO
      DO IUNIT = HIOPEN+NREOQ+1, HIOPEN+NREOQ*(NPROC+1)
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE ( UNIT=IUNIT, STATUS='delete' )
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLTAB ( RTYPE, OQI, IVTYP, MIP, IRQ, BLKND,           40.51
     &                      XC   , YC , XP   , YP )                       40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!
!  2. Purpose
!
!     Printing of table output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     IRQ         request number
!     IVTYP       type of variable output
!     MIP         total number of output points
!     OQI         array containing output request data
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     XP          user x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!     YP          user y-coordinate of output point
!
      INTEGER   IRQ, MIP, OQI(4), IVTYP(OQI(3))
      REAL      BLKND(MXCGL,MYCGL), XC(MIP), YC(MIP), XP(MIP), YP(MIP)    40.51
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EXIST :     logical whether a file exist or not
!     FSTR  :     an auxiliary string
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     LFIELD:     actual length of a part of field OUTLIN
!     MSGSTR:     string to pass message to call MSGERR
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     NUMDEC:     number of decimals in the table
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     OUTLIN:     output line
!
      INTEGER       IENT, IF, IL, ILPOS, IP, IPROC, IUNIT,
     &              IUT, IVTYPE, IXK, IYK, JVAR, LFIELD, NLINES,
     &              NREF, NUMDEC, NVAR
      LOGICAL       EXIST, OPENED
      CHARACTER*80  MSGSTR
      CHARACTER*18  FSTR
      CHARACTER*512 OUTLIN
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read output data from the proper process file and write             40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLTAB')

      NREF   = OQI(1)
      NVAR   = OQI(3)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of heading

            IF ( RTYPE.NE.'TABD' ) THEN
               IF ( RTYPE.EQ.'TABP' .OR. RTYPE.EQ.'TABI' ) THEN
                  NLINES = NLINES + 7
               ELSE IF ( RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABS' ) THEN
                  NLINES = NLINES + 5
                  IF ( RTYPE.EQ.'TABT' ) THEN
                     NLINES = NLINES + 1
                  ELSE
                     IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2
                     NLINES = NLINES + 2 + MIP
                  END IF
                  DO JVAR = 1, NVAR
                     IVTYPE = IVTYP(JVAR)
                     IF ( OVSVTY(IVTYPE).LE.2 ) THEN
                        NLINES = NLINES + 3
                     ELSE
                        NLINES = NLINES + 6
                    END IF
                  END DO
               END IF
            END IF

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

!     --- read output data from the proper process file and write         40.51
!         it in appropriate manner to generic output file                 40.51

      IF ( NREF.NE.PRINTF ) THEN
         IF ( RTYPE.EQ.'TABS' .AND. NSTATM.EQ.1 ) THEN
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               READ (IUNIT,'(A)') OUTLIN
               CALL TXPBLA(OUTLIN,IF,IL)
               IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO
         END IF
         IPLOOP : DO IP = 1, MIP                                          40.51
            IXK = NINT(XC(IP)+100.)-99                                    41.07 40.51
            IYK = NINT(YC(IP)+100.)-99                                    41.07 40.51
            IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.            41.07 40.51
     &           IYK.GT.MYCGL ) THEN                                      41.07 40.51
               CALL WREXCV                                                40.51
               CYCLE IPLOOP                                               40.51
            END IF                                                        40.51
            IPROC = 1                                                     40.51
            PROCLOOP : DO                                                 40.51
              IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                   40.51
                 IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC         40.51
                 READ (IUNIT,'(A)') OUTLIN                                40.51
                 CALL TXPBLA(OUTLIN,IF,IL)                                40.51
                 WRITE (NREF, '(A)') OUTLIN(1:IL)                         40.51
                 EXIT PROCLOOP                                            40.51
              ELSE                                                        40.51
                 IPROC = IPROC + 1                                        40.51
                 IF ( IPROC.LE.NPROC ) THEN                               40.51
                    CYCLE PROCLOOP                                        40.51
                 ELSE                                                     40.51
                    CALL WREXCV                                           40.51
                    EXIT PROCLOOP                                         40.51
                 END IF                                                   40.51
              END IF                                                      40.51
            END DO PROCLOOP                                               40.51
         END DO IPLOOP                                                    40.51
      END IF                                                              40.51
      RETURN

      CONTAINS                                                            40.51
      SUBROUTINE WREXCV                                                   40.51
      IL = 1                                                              40.51
      OUTLIN = '    '                                                     40.51
      IF (RTYPE.EQ.'TABI') THEN                                           40.51
!        --- write point sequence number as first column                  40.51
         WRITE (OUTLIN(1:8), '(I8)') IP                                   40.51
         IL = 9                                                           40.51
         OUTLIN(IL:IL) = ' '                                              40.51
      END IF                                                              40.51
      DO JVAR = 1, NVAR                                                   40.51
         IVTYPE = IVTYP(JVAR)                                             40.51
         IF (IVTYPE.EQ.40) THEN                                           40.51
!           --- for time, 18 characters are needed                        40.51
            FSTR   = '(A18)'                                              40.51
            LFIELD = 18                                                   40.51
            OUTLIN(IL:IL+LFIELD-1) = CHTIME                               40.51
         ELSE                                                             40.51
           IF (RTYPE.EQ.'TABD') THEN                                      40.51
              FSTR   = FLT_TABLE                                          40.51
              LFIELD = FLD_TABLE                                          40.51
           ELSE                                                           40.51
              FSTR   = '(F11.X)'                                          40.51
              LFIELD = 11                                                 40.51
              NUMDEC = MAX (0,6-NINT(LOG10(ABS(OVHEXP(IVTYPE)))))         40.51
              IF (NUMDEC.GT.9) NUMDEC = 9                                 40.51
              WRITE (FSTR(6:6), '(I1)') NUMDEC                            40.51
           END IF                                                         40.51
!          --- write value into OUTLIN                                    40.51
           IF (IVTYPE.EQ.1) THEN                                          40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) XP(IP)             40.51
           ELSE IF (IVTYPE.EQ.2) THEN                                     40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) YP(IP)             40.51
           ELSE                                                           40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) OVEXCV(IVTYPE)     40.51
           END IF                                                         40.51
           IF (OVSVTY(IVTYPE).EQ.3) THEN                                  40.51
              IL = IL + LFIELD + 1                                        40.51
!             --- write second component of a vectorial quantity          40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) OVEXCV(IVTYPE)     40.51
           END IF                                                         40.51
         END IF                                                           40.51
         IL = IL + LFIELD + 1                                             40.51
         OUTLIN(IL-1:IL) = '  '                                           40.51
      END DO                                                              40.51
      CALL TXPBLA(OUTLIN,IF,IL)                                           40.51
      WRITE (NREF, '(A)') OUTLIN(1:IL)                                    40.51
      RETURN                                                              40.51
      END SUBROUTINE WREXCV                                               40.51

      END
!****************************************************************
!
      SUBROUTINE SWCOLSPC ( RTYPE, OQI, MIP, IRQ, BLKND, XC, YC )         40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!
!  2. Purpose
!
!     Printing of spectral output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     IRQ         request number
!     MIP         total number of output points
!     OQI         array containing output request data
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!
      INTEGER   IRQ, MIP, OQI(4)
      REAL      BLKND(MXCGL,MYCGL), XC(MIP), YC(MIP)                      40.51
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EMPTY :     logical whether a line is empty or not
!     EXIST :     logical whether a file exist or not
!     IBLKN :     integer giving node number per subdomain
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IS    :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IXK   :     loop counter
!     IYK   :     loop counter
!     MSGSTR:     string to pass message to call MSGERR
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     OPENED:     logical whether a file is open or not
!     OTYPE :     integer indicating dimension of spectrum
!     OUTLIN:     output line
!
      INTEGER       IBLKN, IENT, IF, IL, ILPOS, IP, IPROC, IS,
     &              IUNIT, IUT, IXK, IYK, NLINES, NREF, OTYPE
      LOGICAL       EMPTY, EXIST, OPENED
      CHARACTER*80  MSGSTR
      CHARACTER (LEN=LENSPO) OUTLIN
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read output data from the proper process file and write             40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLSPC')

      NREF   = OQI(1)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of first part of heading

            NLINES = NLINES + 5
            IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

!           --- write coordinates of output points to generic
!               output file

            DO IP = 1, MIP                                                40.51
               IF ( XC(IP).LT.-0.01 .OR. YC(IP).LT.-0.01 .OR.             40.51
     &              XC(IP).GT.REAL(MXCGL-1)+0.01 .OR.                     40.51
     &              YC(IP).GT.REAL(MYCGL-1)+0.01 ) THEN                   40.51
                  IBLKN = NPROC+1                                         40.51
               ELSE                                                       40.51
                  IXK   = NINT(XC(IP)+100.)-99                            41.07 40.51
                  IYK   = NINT(YC(IP)+100.)-99                            41.07 40.51
                  IBLKN = NINT(BLKND(IXK,IYK))                            40.51
               END IF                                                     40.51
               EMPTY = .TRUE.
               DO IPROC = 1, NPROC
                  IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (EMPTY.AND.IBLKN.EQ.IPROC) THEN                      40.51
                     WRITE (NREF, '(A)') OUTLIN(1:IL)
                     EMPTY = .FALSE.
                  END IF
               END DO
               IF (EMPTY) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO

!           --- count lines of rest of heading and write heading
!               to generic output file

            NLINES = 2 + MSC
            IF (RTYPE(4:4).EQ.'C') THEN
               NLINES = NLINES + 7 + MDC
            ELSE
               NLINES = NLINES + 11
            END IF
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

      IF (RTYPE(4:4).EQ.'C') THEN
         IF (RTYPE.EQ.'SPEC') THEN
            OTYPE = -2
         ELSE
            OTYPE =  2
         END IF
      ELSE
         IF (RTYPE.EQ.'SPE1') THEN
            OTYPE = -1
         ELSE
            OTYPE =  1
         END IF
      END IF

!     --- read output data from the proper process file and write         40.51
!         it in appropriate manner to generic output file                 40.51

      IF ( NSTATM.EQ.1 ) THEN
         DO IPROC = 1, NPROC
            IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
            READ (IUNIT,'(A)') OUTLIN
            CALL TXPBLA(OUTLIN,IF,IL)
            IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
         END DO
      END IF

      IPLOOP : DO IP = 1, MIP                                             40.51
         IXK = NINT(XC(IP)+100.)-99                                       41.07 40.51
         IYK = NINT(YC(IP)+100.)-99                                       41.07 40.51
         IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.               41.07 40.51
     &        IYK.GT.MYCGL ) THEN                                         41.07 40.51
            WRITE (NREF, '(A)') 'NODATA'                                  40.51
            CYCLE IPLOOP                                                  40.51
         END IF                                                           40.51
         IPROC = 1                                                        40.51
         PROCLOOP : DO                                                    40.51
           IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                      40.51
              IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC            40.51
              READ (IUNIT,'(A)') OUTLIN                                   40.51
              CALL TXPBLA(OUTLIN,IF,IL)                                   40.51
              WRITE (NREF, '(A)') OUTLIN(1:IL)                            40.51
              IF ( OUTLIN(IF:IL).NE.'NODATA' ) THEN                       40.51
                 IF ( ABS(OTYPE).EQ.1 ) THEN                              40.51
                    DO IS = 1, MSC                                        40.51
                       READ (IUNIT,'(A)') OUTLIN                          40.51
                       CALL TXPBLA(OUTLIN,IF,IL)                          40.51
                       WRITE (NREF, '(A)') OUTLIN(1:IL)                   40.51
                    END DO                                                40.51
                 ELSE                                                     40.51
                    IF ( OUTLIN(IF:IL).NE.'ZERO' ) THEN                   40.51
                       DO IS = 1, 1+MSC                                   40.51
                          READ (IUNIT,'(A)') OUTLIN                       40.51
                          CALL TXPBLA(OUTLIN,IF,IL)                       40.51
                          WRITE (NREF, '(A)') OUTLIN(1:IL)                40.51
                       END DO                                             40.51
                    END IF                                                40.51
                 END IF                                                   40.51
              END IF                                                      40.51
              EXIT PROCLOOP                                               40.51
           ELSE                                                           40.51
              IPROC = IPROC + 1                                           40.51
              IF ( IPROC.LE.NPROC ) THEN                                  40.51
                 CYCLE PROCLOOP                                           40.51
              ELSE                                                        40.51
                 WRITE (NREF, '(A)') 'NODATA'                             40.51
                 EXIT PROCLOOP                                            40.51
              END IF                                                      40.51
           END IF                                                         40.51
         END DO PROCLOOP                                                  40.51
      END DO IPLOOP

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLBLK ( RTYPE, OQI, IVTYP, FAC  , PSNAME,
     &                      MXK  , MYK, IRQ  , BLKND, XC    ,             40.51
     &                      YC   )                                        40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE OUTP_DATA
      USE M_PARALL
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
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
!
!  0. Authors
!
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.31, Dec. 03: New subroutine
!     40.41, Jun. 04: some improvements with respect to MATLAB
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!
!  2. Purpose
!
!     Writing of block output by means of collecting
!     individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     FAC         factors of multiplication of block output
!     IRQ         request number
!     IVTYP       type of variable output
!     MXK         number of points in x-direction of output frame
!     MYK         number of points in y-direction of output frame
!     OQI         array containing output request data
!     PSNAME      name of output locations
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!
      INTEGER   MXK, MYK, IRQ, OQI(4), IVTYP(OQI(3))
      REAL      BLKND(MXCGL,MYCGL), XC(MXK*MYK), YC(MXK*MYK)              40.51
      REAL      FAC(OQI(3))
      CHARACTER RTYPE*4, PSNAME*8
!
!  6. Local variables
!
!     CTIM  :     string representing date of computation
!     DFAC  :     multiplication factor of block output
!     EXIST :     logical whether a file exist or not
!     FMAX  :     auxiliary real
!     FTIP  :     auxiliary real
!     FTIP1 :     auxiliary real
!     FTIP2 :     auxiliary real
!     IDLA  :     lay-out indicator
!     IF    :     first non-character in string
!     IFAC  :     auxiliary integer
!     IENT  :     number of entries
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPD   :     switch for printing on paper or writing to file
!     IPROC :     loop counter
!     IREC  :     direct access file record counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     MATLAB:     indicates whether binary Matlab files are used
!     MSGSTR:     string to pass message to call MSGERR
!     NAMVAR:     name of MATLAB variable
!     NREF  :     unit reference number
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     VOQ   :     collected output variables
!
      INTEGER      IDLA, IENT, IF, IFAC, IL, ILPOS, IP, IPD,
     &             IPROC, IUNIT, IUT, IXK, IYK, IVTYPE, JVAR, NREF,
     &             NVAR
      REAL         DFAC, FMAX, FTIP, FTIP1, FTIP2
      LOGICAL      EXIST, OPENED
      INTEGER, SAVE :: IREC(MAX_OUTP_REQ)=0                               40.51 40.41
      LOGICAL, SAVE :: MATLAB=.FALSE.                                     40.41
      CHARACTER*80 MSGSTR
      CHARACTER (LEN=20) :: CTIM                                          40.41
      CHARACTER (LEN=30) :: NAMVAR                                        40.41
      REAL, ALLOCATABLE :: VOQ(:,:)
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     SBLKPT           Writes block output to an ASCII file
!     STRACE           Tracing routine for debugging
!     SWRMAT           Writes block output to a binary Matlab file
!     TABHED           Prints heading
!     TXPBLA           Removes leading and trailing blanks in string      40.41
!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           open individual process output files
!
!     read data from the proper process file and write                    40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLBLK')

      NREF = OQI(1)
      NVAR = OQI(3)
      IDLA = OQI(4)

      IF ( RTYPE.EQ.'BLKP' ) THEN
        IPD = 1
        IF (NREF.EQ.PRINTF) CALL TABHED ('SWAN', PRINTF)
      ELSE IF ( RTYPE.EQ.'BLKD' ) THEN
        IPD = 2
      ELSE
        IPD = 3
      ENDIF

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )
         MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                       40.41
     &            INDEX (FILENM, '.mat' ).NE.0                            40.41

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

            IF (MATLAB .AND. .NOT.OPENED) THEN
               CLOSE(NREF)
               OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',
     &              ACCESS='DIRECT', RECL=4)                              41.08
               IREC(IRQ) = 1                                              40.51
            END IF

!           --- open individual process output files

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
            END DO

         END IF

      END IF

!     --- read data from the proper process file and write                40.51
!         it in appropriate manner to generic output file                 40.51

      CTIM = CHTIME                                                       40.41
      CALL TXPBLA(CTIM,IF,IL)                                             40.41
      CTIM(9:9)='_'                                                       40.41

      ALLOCATE(VOQ(MXK*MYK,2))

      DO JVAR = 1, NVAR

         IVTYPE = IVTYP(JVAR)
         DFAC   = FAC(JVAR)

         VOQ = OVEXCV(IVTYPE)

         DO IP = 1, MXK*MYK                                               40.51
            IXK = NINT(XC(IP)+100.)-99                                    41.07 40.51
            IYK = NINT(YC(IP)+100.)-99                                    41.07 40.51
            IPROC = 1                                                     40.51
            PROCLOOP1 : DO                                                40.51
              IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                   40.51
                 IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC         40.51
                 READ (IUNIT, FLT_BLKP) VOQ(IP,1)                         40.51
                 EXIT PROCLOOP1                                           40.51
              ELSE                                                        40.51
                 IPROC = IPROC + 1                                        40.51
                 IF ( IPROC.LE.NPROC ) THEN                               40.51
                    CYCLE PROCLOOP1                                       40.51
                 ELSE                                                     40.51
                    EXIT PROCLOOP1                                        40.51
                 END IF                                                   40.51
              END IF                                                      40.51
            END DO PROCLOOP1                                              40.51
         END DO                                                           40.51

         IF ( OVSVTY(IVTYPE).GE.3 ) THEN

            DO IP = 1, MXK*MYK                                            40.51
               IXK = NINT(XC(IP)+100.)-99                                 41.07 40.51
               IYK = NINT(YC(IP)+100.)-99                                 41.07 40.51
               IPROC = 1                                                  40.51
               PROCLOOP2 : DO                                             40.51
                 IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                40.51
                    IUNIT=HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC        40.51
                    READ (IUNIT, FLT_BLKP) VOQ(IP,2)                      40.51
                    EXIT PROCLOOP2                                        40.51
                 ELSE                                                     40.51
                    IPROC = IPROC + 1                                     40.51
                    IF ( IPROC.LE.NPROC ) THEN                            40.51
                       CYCLE PROCLOOP2                                    40.51
                    ELSE                                                  40.51
                       EXIT PROCLOOP2                                     40.51
                    END IF                                                40.51
                 END IF                                                   40.51
               END DO PROCLOOP2                                           40.51
            END DO                                                        40.51
         END IF

         IF ( IPD.EQ.1 ) THEN
            IF ( DFAC.LE.0. ) THEN
               IF ( OVHEXP(IVTYPE).LT.0.5E10 ) THEN
                  IFAC = INT (10.+LOG10(OVHEXP(IVTYPE))) - 13
               ELSE
                  IF ( OVSVTY(IVTYPE).EQ.1 ) THEN
                     FMAX = 1.E-8
                     DO IP = 1, MXK*MYK
                        FTIP = ABS(VOQ(IP,1))
                        FMAX = MAX (FMAX, FTIP)
                     END DO
                  ELSE IF ( OVSVTY(IVTYPE).EQ.2 ) THEN
                     FMAX = 1000.
                  ELSE IF ( OVSVTY(IVTYPE).EQ.3 ) THEN
                     FMAX = 1.E-8
                     DO IP = 1, MXK*MYK
                        FTIP1 = ABS(VOQ(IP,1))
                        FTIP2 = ABS(VOQ(IP,2))
                        FMAX  = MAX (FMAX, FTIP1, FTIP2)
                     END DO
                  END IF
                  IFAC = INT (10.+LOG10(FMAX)) - 13
               END IF
               DFAC = 10.**IFAC
            END IF
         ELSE
           IF ( DFAC.LE.0. ) DFAC = 1.
         END IF

         IF (OVSVTY(IVTYPE) .LT. 3) THEN
            IF (MATLAB) THEN
               IF (IL.EQ.1 .OR. IVTYPE.LT.3 .OR. IVTYPE.EQ.52) THEN       40.94 40.41
                  NAMVAR = OVSNAM(IVTYPE)                                 40.41
               ELSE                                                       40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_'//CTIM                                      40.41
               END IF                                                     40.41
               CALL SWRMAT( MYK, MXK, NAMVAR, VOQ(1,1), NREF,
     &                      IREC(IRQ), IDLA, OVEXCV(IVTYPE) )             40.51
            ELSE
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE), VOQ(1,1) )
            END IF
         ELSE
            IF (MATLAB) THEN
               IF (IL.EQ.1) THEN                                          40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_x'                                           40.41
               ELSE                                                       40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_x_'//CTIM                                    40.41
               END IF                                                     40.41
               CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,1), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )  40.51
               IF (IL.EQ.1) THEN                                          40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_y'                                           40.41
               ELSE                                                       40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_y_'//CTIM                                    40.41
               END IF                                                     40.41
               CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,2), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )  40.51
            ELSE
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE)//'X-comp',
     &                      VOQ(1,1) )
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE)//'Y-comp',
     &                      VOQ(1,2) )
            END IF
         END IF

      END DO

      IF (IPD.EQ.1 .AND. NREF.EQ.PRINTF) WRITE (PRINTF, 111)

      DEALLOCATE(VOQ)

  111 FORMAT (///)

      RETURN
      END
