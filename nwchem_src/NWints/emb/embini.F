C This subroutine initialize the embedding potential
C Loading box size/shapes, potential dimensions etc.
      subroutine embini

      implicit none

#include "embcom.fh"
#include "global.fh"
#include "stdio.fh"

C       Input file
        character*6 finput
        parameter (finput='emb.in')

C       I/O stuff
        integer iunit,idx,iost
        character*80 cline

C       local variables
        integer iproc,nproc,i,j
        double precision vec(3)

        lembed = .true.
        iproc = ga_nodeid()
        nproc = ga_nnodes()
        iunit=240
        open(iunit,file=finput,status='old',iostat=iost)
        if ( iost/=0 ) then
          if ( iproc == 0 ) then
            write(luout,*) ' error on opening emb.in, stop!'
          endif
          close(iunit)
          stop
        endif

        if ( iproc == 0 ) then
          write(luout,*) 'Loading parameters from emb.in...'
        endif

C       read grid dimension
        call embget(iunit,'grid_dim',cline,idx)
        if (idx.eq.0) then
          if ( iproc == 0 ) then
            write(luout,*) 'error: cannot read grid_dim in emb.in'
            write(luout,*) 'Please provide:'
            write(luout,*) 'grid_dim: n1 n2 n3'
            write(luout,*) 'exiting...'
          endif
          stop
        endif
        read(iunit,*) n1,n2,n3
        nunipt = n1*n2*n3

C Real-space lattice vectors in Angstrom
        call embget(iunit, 'rlatv', cline, idx)
        if (idx.eq.0) then
          if ( iproc == 0 ) then
            write(luout,*) 'error: cannot read raltc in emb.in'
            write(luout,*) 'Please provide:'
            write(luout,*) 'rlatv:'
            write(luout,*) '(1,1)  (1,2)  (1,3)'
            write(luout,*) '(2,1)  (2,2)  (2,3)'
            write(luout,*) '(3,1)  (3,2)  (3,3)'
            write(luout,*) 'exiting...'
          endif
          stop
        endif
        do i=1,3
          read(iunit,*) (rlatv(i,j),j=1,3)
        end do
        rlatv = rlatv/ANGSTR

CJin: check tag shift:
        call embget(iunit, 'shift', cline, idx)
        if (idx.eq.0) then
          if ( iproc == 0 ) then
            write(luout,*) 'no shifted vector is provided'
          endif
          shift(:) = 0.0
        else
          if ( iproc == 0 ) then
            write(luout,*) 'shifted vector is provded'
          endif
          read(iunit,*) shift
        end if
        shift = shift/ANGSTR

        close(iunit)

C Calculate real-space unit cell volume
        call cross_product( rlatv(1,:), rlatv(2,:), vec(:) )
        vol = dot_product( rlatv(3,:), vec(:) )
        if ( iproc == 0 ) then
          write(6,*) ' rlatv='
          do i=1,3
             write(luout,'(3(TR3,E14.6))') (rlatv(i,j),j=1,3)
          end do
          write(luout,*) ' vol= ', vol
        endif
        volelm = vol / nunipt

      end subroutine embini

C     This subroutine initialize the fully nonlocal embedding
C     calculation, work on density matrix elements space
      subroutine nlembini(rtdb)
      implicit none
#include "embcom.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb, geom, basis, nbf

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

c     set embedding info tag
      nlembed = .true.

c     load geometry and basis set info, need it to create potential
c     matrix
      if (.not. geom_create(geom, 'geometry'))
     $    call errquit('nlembini: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $    call errquit('nlembini: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('nlembini: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('nlembini: no ao basis ', 0, RTDB_ERR)
c     get number of AO basis functions
      if (.not. bas_numbf(basis,nbf)) call
     &      errquit('nlembini: could not get nbf',0, BASIS_ERR)

c     create handle for nonlocal potential matrix
      g_pot = ga_create_atom_blocked( geom, basis, 'nonloc Vemb' )

c     read in the nonlocal embedding potential
      call ao_1prdm_read_asicii(nbf, g_pot, 'nlembpot.dat')

      if (ga_nodeid().eq.0) then
        write(luout,*) ' Nonlocal Vemb matrix loaded successfully.'
      endif 

      return 

      end subroutine nlembini


      subroutine embget(iunit,key,cline,idx)

      implicit none

C     ------------------------------------------------------------------
C     Inputs
C     Input file unit
      integer iunit
C     Key to search for
      character*(*) key

C     Outputs
C     Line containing key and index to the first position after key,
C     idx=0 on output if key is not found
      character*80 cline
      integer idx
C     ------------------------------------------------------------------

      integer ios

C     Make sure we are at the top of the file
      rewind iunit

      idx = 0
      ios = 0
      do while ((idx.eq.0).and.(ios.eq.0))
C     Clear character line buffer and read a line
         cline = ''
         read(iunit,'(A80)',iostat=ios) cline

C     Look for occurance of key, ignore lines beginning with #
         if (cline(1:1).ne.'#') then
            idx = index(cline,key)
         end if
      end do

      if (idx.gt.0) then
         idx = idx+len(key)
      end if

      end subroutine embget

      subroutine embunset
        implicit none
#include "embcom.fh"
        lembed = .false.
        nlembed = .false.

      end subroutine embunset
