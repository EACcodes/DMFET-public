c
c     Write out 1-particle reduced density matrix
c
      subroutine ao_1prdm_write(nbf,g_dens,filename)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer nbf         ! [input] Number of basis functions
      integer g_dens      ! [input] 1-particle reduced density matrix
      character*255 filename
c
      integer l_dens, k_dens
      integer ok
      integer inntsize
      integer i
      character*32 pname
c
      integer unitno
      parameter (unitno = 78)
c
c     Preliminaries
      pname = "ao_1prdm_write: "
      l_dens = -1               ! An invalid MA handle
c
c     Read the file
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
       open(unitno,status='unknown',form='unformatted',file=filename)
       if (.not. ma_push_get(mt_dbl,nbf,'ao_1prdm_write',
     &     l_dens,k_dens)) call errquit(pname//'ma failed', nbf, MA_ERR)
       do i = 1,nbf
            call ga_get(g_dens, 1, nbf, i, i, dbl_mb(k_dens),1)
            call swrite(unitno, dbl_mb(k_dens), nbf)
       enddo
       close(unitno)
c
c      Clean up
       if (.not. ma_pop_stack(l_dens))
     $      call errquit(pname//'ma pop failed', l_dens, MA_ERR)
      endif   ! ga_nodeid()
c
c     Broadcast status to other nodes
      ok = 1
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
      call ga_sync()
c
      end
c $Id: ao_1prdm_write.F 26040 2014-08-26 03:43:05Z niri $

c--------------------------------------------------------------------------------
c
c     Write out 1-particle reduced density matrix, in asic ii form
c     added by kuang, modified from ao_1pdrm_write
c
      subroutine ao_1prdm_write_asicii(nbf,g_dens,filename)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer nbf         ! [input] Number of basis functions
      integer g_dens      ! [input] 1-particle reduced density matrix
      character*255 filename
c
      integer l_dens, k_dens
      integer ok
      integer inntsize
      integer i, j
      character*32 pname
c
      integer unitno
      parameter (unitno = 78)
c
c     Preliminaries
      pname = "ao_1prdm_write: "
      l_dens = -1               ! An invalid MA handle
c
c     Read the file
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
       open(unitno,status='unknown',file=filename)
       if (.not. ma_push_get(mt_dbl,nbf,'ao_1prdm_write',
     &     l_dens,k_dens)) call errquit(pname//'ma failed', nbf, MA_ERR)
       do i = 1,nbf
            call ga_get(g_dens, 1, nbf, i, i, dbl_mb(k_dens),1)
            do j = 1,nbf
               write(unitno,'(E25.15)',advance='no') dbl_mb(k_dens+j-1)
               if ( mod(j,5).eq.0 ) write(unitno,'(A)')
     &                              '' ! change line
            enddo
            !if ( mod(j-1,5).ne.0 ) write(unitno,'(A)') '' ! change line
            write(unitno,'(A)') '' ! change line
       enddo
       close(unitno)
c
c      Clean up
       if (.not. ma_pop_stack(l_dens))
     $      call errquit(pname//'ma pop failed', l_dens, MA_ERR)
      endif   ! ga_nodeid()
c
c     Broadcast status to other nodes
      ok = 1
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
      call ga_sync()
c
      end
c $Id: ao_1prdm_write.F 26040 2014-08-26 03:43:05Z niri $
