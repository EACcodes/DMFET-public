C============================================================
C This subroutine computes the integral due to embedding
C potential.
C============================================================

      Subroutine embint(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &                  Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &                  Nints,NPP,bNAI,canAB)
        implicit none
#include "embcom.fh"
#include "global.fh"
#include "stdio.fh"

C       input variables
        logical canAB
        integer NPA,NPB,NCA,NCB,La,Lb,Nints,NPP
        double precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
        double precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)
        double precision bNAI(Nints)

C       local variabels
        integer nn, La2, Lb2, iproc, ioflag
        integer a, b, b_limit, m, n
        logical :: flag_read, flag_fileexist, io_node
C       primitive integration matrix
        double precision Final(NPA,NPB,(La+1)*(La+2)/2,(Lb+1)*(Lb+2)/2)
C       tmp variables for io
        character(len=300) stmp
        integer tmp_La,tmp_Lb,i,iost
        double precision tmp_Axyz(3),tmp_Bxyz(3), tmp_Aexp, tmp_Bexp
        double precision :: tmp
        character(100) :: buffer(4)

        if ( .not. lembed ) then
          return
        endif

        La2 = (La+1)*(La+2)/2
        Lb2 = (Lb+1)*(Lb+2)/2

        iproc = ga_nodeid()
        io_node = .false.
        if (iproc==0) io_node = .true.

        flag_read = .false.
        ! Whether these shells have be computed or not
        ! Check the existing oneint.dat file
        inquire(FILE="oneint.dat", exist=flag_fileexist)
        ! the reading would only work in the serial run.
        ! in the parallel run, the format of oneint.dat will be messed
        ! up.
        if ( flag_fileexist .and. ga_nnodes()==1 ) then
          open(unit=111,file='oneint.dat',action='read')
          ioflag = 0
          do while (ioflag >= 0)
            !write(luout,*) 'Read oneint.dat.'
            read(111,'(A13,I3,A4,I3,E13.4,E13.4,A5,3F7.3,A7,3F7.3,A3,F)'
     $          ,iostat=ioflag) buffer(1),tmp_La,buffer(2),tmp_Lb
     $          ,tmp_Aexp, tmp_Bexp
     $          ,buffer(3),tmp_Axyz(1),tmp_Axyz(2),tmp_Axyz(3)
     $          ,buffer(4),tmp_Bxyz(1),tmp_Bxyz(2),tmp_Bxyz(3)
     $          ,buffer(5),tmp
            ! found the right block in the oneint.dat file
            if(tmp_La==La .and. tmp_Lb==Lb .and.
     $         abs(Aprims(1)-tmp_Aexp)/Aprims(1) < 0.01 .and.
     $         abs(Bprims(1)-tmp_Bexp)/Bprims(1) < 0.01 .and.
     $         dot_product(tmp_Axyz-Axyz,tmp_Axyz-Axyz)<1d-4 .and.
     $         dot_product(tmp_Bxyz-Bxyz,tmp_Bxyz-Bxyz)<1d-4 ) then
              flag_read = .true.
              ! start reading the block
              do b = 1, Lb2
                do a = 1, La2
                  do n = 1, NPB
                    do m = 1, NPA
                      Final(m,n,a,b) = tmp
                      read(111,'(A13,I3,A4,I3,A31,3F7.3,A7,3F7.3,A3,F)'
     $                 ,iostat=ioflag)
     $                 ,buffer(1),tmp_La,buffer(2),tmp_Lb
     $                 ,buffer(3),tmp_Axyz(1),tmp_Axyz(2),tmp_Axyz(3)
     $                 ,buffer(4),tmp_Bxyz(1),tmp_Bxyz(2),tmp_Bxyz(3)
     $                 ,buffer(5),tmp
                    enddo
                  enddo
                enddo
              enddo ! end of reading the block
              goto 222
            endif
          enddo
222       continue
          close(111)
        endif
      
        ! file doesn't exist or file does not include the right block
        ! then do the real calculation.
        if(flag_fileexist==.false. .or. flag_read==.false.) then
          call uniint(Aprims,NPA,Bprims,NPB,Axyz,Bxyz,La,Lb,Final)
        endif
        ! now we have primitive integrals, write them to oneint.dat
        ! again, the file i/o only works with serial run
        if(flag_read==.false. .and. ga_nnodes()==1) then
          open(unit=111,file='oneint.dat',
     &       action='write',access='append')
          do b=1,Lb2
            do a=1,La2
              do n=1,NPB
                do m=1,NPA
                  if (ABS(Final(m,n,a,b))<1e-30) Final(m,n,a,b) = 0.d0
                write(111,
     $          '(a,a,i3,a,i3,e12.4,e12.4,2(a,3(f7.3)),a,es22.12)')
     $       ' embint2:',' la=',la,' lb=',lb,
     $       Aprims(m),Bprims(n),
     $       ' Apos=',Axyz(1),Axyz(2),Axyz(3),
     $       ' RBpos=',Bxyz(1),Bxyz(2),Bxyz(3),' => ',Final(m,n,a,b)
                end do
              end do
            end do
          end do
          close(111)
        endif

        ! now use the primitive integrals to compute the contraction
        ! integrals
        nn = 0
        do a = 1,La2
          if( canAB )then
            b_limit = a
          else
            b_limit = Lb2
          end if
          do b = 1, b_limit
            tmp = 0.0d0
            do m = 1, NPA
              do n = 1, NPB
                ! note, no gencon is supported, NCA=NCB=1 by default
                ! multiply by contraction coefficients
                tmp = tmp + Final(m,n,a,b)*Acoefs(m,1)*Bcoefs(n,1)
              enddo 
            enddo
            nn = nn + 1
            ! add the one electron integrals to the <a|V|b> term
            bNAI(nn) = bNAI(nn) + tmp
          enddo
        enddo

        return
 
      end subroutine

C============================================================
C Integrate with respect to the embedding potential potuni and add
C contribution to array Final
C============================================================
      subroutine uniint(expA,nexpA,expB,nexpB,Apos,RBpos,la,lb,Final)
        implicit none
#include "embcom.fh"
#include "global.fh"
#include "stdio.fh"
#include "fftw3.f"
#include "mafdecls.fh"
#include "errquit.fh"
        integer nexpA,nexpB,la,lb
C       double precision expA(nexpA),expB(nexpB),Apos(3),RBpos(3),
C    $       Final(nexpA,nexpB,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
        double precision expA(nexpA),expB(nexpB),Apos(3),RBpos(3),
     $       Final(nexpA,nexpB,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)

C       Work arrays
        integer mexp,nexp

        double precision x,y,z,xA,yA,zA,xRB,yRB,zRB,rAsq,rRBsq,foo
        double precision center(3)
        integer ncartA,ncartB
        integer i,m,n,a,b,ix,iy,iz,ioper,noper
CPeilin
        integer IPCOUNT
        integer ndumind,ndumtot
C       integer nembdum
        double precision embnoper
        double precision dumprint

        double precision tmp
CPeilin
CKuang
        integer nproc,iproc
        integer ipt,i1,j1,k1,i2,j2,k2,k3,imax,jmax,kmax,npts
        integer ixlo,ixhi,iylo,iyhi,izlo,izhi,ixcell,iycell,izcell
        integer mlo, alo
        integer ndumind1
        integer*8 :: fft_plan

        integer :: ivq, lvq, ivr, lvr, ikx, lkx, iky, lky, ikz, lkz
        integer :: iphase, lphase
        integer :: i1cart, j1cart, k1cart, i2cart, j2cart, k2cart
        integer :: li1cart, lj1cart, lk1cart, li2cart, lj2cart, lk2cart
        integer :: iacart, ibcart, lacart, lbcart

!        double complex, allocatable :: vq(:,:,:), phase(:)
!        double precision, allocatable :: vr(:,:,:), kx(:),ky(:),kz(:)
        double precision :: exp_thresh=20d0, alpha_cut=20.0d0 !a.u.
        double precision :: box(3,3), box_inv(3,3)
!        integer, allocatable ::  i1cart(:),j1cart(:),k1cart(:)
!        integer, allocatable ::  i2cart(:),j2cart(:),k2cart(:)
!        double precision, allocatable ::  acart(:), bcart(:)
        double precision :: fac, kRI, alpha1,alpha2,rrbound
        double precision :: RI(3),r(3),dr(3),ds(3),dsp(3)
        double complex :: f1(0:la+lb), f2(0:la+lb), f3(0:la+lb)
        logical :: lsamecenter, lsharp, lsameshell, io_node
        save ivq, ivr, ikx, iky, ikz, iphase, box, box_inv
        save lvq, lvr, lkx, lky, lkz, lphase

C       Number of basis function in shells A and B
        ncartA = (la+1)*(la+2)/2
        ncartB = (lb+1)*(lb+2)/2
        ndumtot=ncartB*ncartA*nexpB*nexpA

        npts = (n1/2+1)*n2*n3

        ! prepare the fft array
        if ( .not. ma_get_index(lvq, ivq) ) then
!          allocate( vq(n1/2+1,n2,n3), vr(n1, n2, n3) )
!          allocate( phase(npts), kx(npts),ky(npts),kz(npts) )
          if ( .not. ma_alloc_get(mt_dcpl, npts, 'vq', lvq, ivq) )
     &      call errquit('uniint: cannot allocate vq',0, MA_ERR)
          if ( .not. ma_alloc_get(mt_dbl, nunipt, 'vr', lvr, ivr) )
     &      call errquit('uniint: cannot allocate vr',0, MA_ERR)
          if ( .not.ma_alloc_get(mt_dcpl, npts, 'phase', lphase,iphase))
     &      call errquit('uniint: cannot allocate phase',0, MA_ERR)
          if ( .not. ma_alloc_get(mt_dbl, npts, 'kx', lkx, ikx) )
     &      call errquit('uniint: cannot allocate kx',0, MA_ERR)
          if ( .not. ma_alloc_get(mt_dbl, npts, 'ky', lky, iky) )
     &      call errquit('uniint: cannot allocate ky',0, MA_ERR)
          if ( .not. ma_alloc_get(mt_dbl, npts, 'kz', lkz, ikz) )
     &      call errquit('uniint: cannot allocate kz',0, MA_ERR)
          open(unit=111,file='embpot.dat',action='read')
          i = 1
          do iz = 1, n3
            do iy = 1, n2
              do ix = 1, n1
                read(111,*) dbl_mb(ivr+i-1) 
                i = i+1
              enddo
            enddo
          enddo
          close(111)
          ! Fourier transform of external potential
          !vq = 0.0d0
          call dfftw_plan_dft_r2c_3d(fft_plan, n1, n2, n3, dbl_mb(ivr)
     $                               , dcpl_mb(ivq), FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c(fft_plan,dbl_mb(ivr),dcpl_mb(ivq))
          do i = 1, npts
            dcpl_mb(ivq+i-1) = dcpl_mb(ivq+i-1)/nunipt
          enddo
          call dfftw_destroy_plan(fft_plan)
          ! setup box
          box = transpose( rlatv ) ! change to column-based layout
          ! construct reciprocal box
          call cross_product(box(:,1),box(:,2),box_inv(:,3))
          call cross_product(box(:,2),box(:,3),box_inv(:,1))
          call cross_product(box(:,3),box(:,1),box_inv(:,2))
          ! box_inv contains a*,b*,c* on column
          ! actually it is transpose(box_inv)
          box_inv = box_inv/vol
          ! setup k-points
          call setup_kpoints(box, box_inv, n1, n2, n3, npts,
     $                       dbl_mb(ikx), dbl_mb(iky), dbl_mb(ikz))
        endif
C Volume element
        volelm = vol/dble(nunipt)

C       Allocate scratch arrays
!        allocate( i1cart(ncartA), j1cart(ncartA), k1cart(ncartA) )
!        allocate( i2cart(ncartB), j2cart(ncartB), k2cart(ncartB) )
!        allocate( acart(ncartA), bcart(ncartB) )
        if (.not.ma_push_get(mt_int, ncartA, 'i1cart', li1cart, i1cart))
     &    call errquit('uniint: cannot allocate i1cart',0, MA_ERR)
        if (.not.ma_push_get(mt_int, ncartA, 'j1cart', lj1cart, j1cart))
     &    call errquit('uniint: cannot allocate i1cart',0, MA_ERR)
        if (.not.ma_push_get(mt_int, ncartA, 'k1cart', lk1cart, k1cart))
     &    call errquit('uniint: cannot allocate i1cart',0, MA_ERR)
        if (.not.ma_push_get(mt_int, ncartB, 'i2cart', li2cart, i2cart))
     &    call errquit('uniint: cannot allocate i2cart',0, MA_ERR)
        if (.not.ma_push_get(mt_int, ncartB, 'j2cart', lj2cart, j2cart))
     &    call errquit('uniint: cannot allocate j2cart',0, MA_ERR)
        if (.not.ma_push_get(mt_int, ncartB, 'k2cart', lk2cart, k2cart))
     &    call errquit('uniint: cannot allocate k2cart',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl, ncartA, 'acart', lacart, iacart))
     &    call errquit('uniint: cannot allocate acart',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl, ncartB, 'bcart', lbcart, ibcart))
     &    call errquit('uniint: cannot allocate bcart',0, MA_ERR)

        Final = 0.0d0
        ! whether the two shells belongs to the same center?
        if (dot_product(Apos-RBpos,Apos-RBpos)<1.0d-6) then
          lsamecenter = .true.
        else
          lsamecenter = .false.
        endif

        ! whether it is same shell? 
        ! the transpose symmetry disabled in NWchem.
        lsameshell = .false.
        !if (lsamecenter .and. la==lb) then
        !  lsameshell = .true.
        !else
        !  lsameshell = .false.
        !endif

        ! setup the angular index
        call set_angular_index( la, int_mb(i1cart), int_mb(j1cart)
     &                            , int_mb(k1cart), ncartA )
        call set_angular_index( lb, int_mb(i2cart), int_mb(j2cart)
     &                            , int_mb(k2cart), ncartB )
        ! if the two shells are not in the same center, skip reciprocal
        ! part of the code
        if (.not. lsamecenter) goto 888
        ! setup phases
        ipt = 1
        do iz = 1, n3
          do iy = 1, n2
            do ix = 1, n1/2+1
              kRI = dbl_mb(ikx+ipt-1)*Apos(1) +dbl_mb(iky+ipt-1)*Apos(2)
     &            + dbl_mb(ikz+ipt-1)*Apos(3)
              dcpl_mb(iphase+ipt-1) = cmplx(cos(kRI), sin(kRI)) ! e^i*k*RI
              ipt = ipt + 1
            enddo
          enddo
        enddo

        do n = 1, nexpB
          if (lsameshell) then
            mlo = n
          else
            mlo = 1
          endif
          do m = mlo, nexpA
            ! if the integrand is not sharp enough, skip
            if (expA(m)+expB(n)<=alpha_cut) then
              goto 777
            endif
            ipt = 1
            do iz = 1, n3
              do iy = 1, n2
                do ix = 1, n1/2+1
                  !if (nproc>1) then
                  !  ! this point is not this process's job, skip
                  !  if ( mod(ipt,nproc) /= iproc ) goto 222
                  !endif
                  call int_pw_gauss(dbl_mb(ikz+ipt-1), la+lb, 
     &                                     expA(m)+expB(n), f3)
                  call int_pw_gauss(dbl_mb(iky+ipt-1), la+lb, 
     &                                     expA(m)+expB(n), f2)
                  call int_pw_gauss(dbl_mb(ikx+ipt-1), la+lb, 
     &                                     expA(m)+expB(n), f1)
                  if ( ix-1 == 0 .or. n1 - (ix-1) == ix-1) then
                    fac = 1.0
                  else
                    fac = 2.0
                  endif
                  do b = 1, ncartB
                    i2 = int_mb(i2cart+b-1)
                    j2 = int_mb(j2cart+b-1)
                    k2 = int_mb(k2cart+b-1)
                    if (lsameshell .and. m==n) then
                      alo = b
                    else
                      alo = 1
                    endif
                    do a = alo, ncartA
                      i1 = int_mb(i1cart+a-1)
                      j1 = int_mb(j1cart+a-1)
                      k1 = int_mb(k1cart+a-1)
                      Final(m,n,a,b) = Final(m,n,a,b) + fac*
     $                            REAL(dcpl_mb(ivq+ipt-1) 
     $                               * dcpl_mb(iphase+ipt-1)
     $                            *f1(i1+i2)*f2(j1+j2)*f3(k1+k2))
                    enddo
                  enddo
222               continue
                  ipt = ipt + 1
                enddo
              enddo
            enddo
777         continue
          enddo
        enddo

888     continue
        ! Now take care of the real space integration for softer
        ! off-center Gaussian pairs

        do n = 1, nexpB
          alpha2 = expB(n)
          if ( lsameshell ) then
            mlo = n
          else
            mlo = 1
          endif
          do m = mlo, nexpA
            alpha1 = expA(m)
            ! check whether this paire has been computed in the
            ! reciprocal integration or not
            if (lsamecenter==.true. .and. alpha1+alpha2>alpha_cut) then
              goto 999 ! skip pair
            endif
           fac = dot_product(Apos-RBpos, Apos-RBpos)
     $            * alpha1*alpha2 / (alpha1+alpha2)
            ! check the size of the prefactor
            if ( fac > exp_thresh ) then !skip pair
              goto 999 ! skip pair
            else
              fac = exp(-fac)
            endif
            ! Center of the combined Gaussian
            RI = (Apos(:)*alpha1+RBpos(:)*alpha2) / (alpha1+alpha2)
            ! A little fail-safe procedure for the diffuse basis
            ! Some diffuse basis may not be contained in a single cell
            rrbound = 0.25d0/dot_product(box_inv(:,1),box_inv(:,1))
            call set_ibound( rrbound, alpha1+alpha2, ixlo, ixhi
     $                       , exp_thresh )
            rrbound = 0.25d0/dot_product(box_inv(:,2),box_inv(:,2))
            call set_ibound( rrbound, alpha1+alpha2, iylo, iyhi
     $                       , exp_thresh )
            rrbound = 0.25d0/dot_product(box_inv(:,3),box_inv(:,3))
            call set_ibound( rrbound, alpha1+alpha2, izlo, izhi
     $                       , exp_thresh )
            ix = 0
            iy = 0
            iz = 0
            do i = 1, nunipt
              !if ( nproc > 1 ) then
              !  ! this point is not this process's job,
              !  ! skip
              !  if (mod(i,nproc)/=iproc) goto 333
              !endif
              r(1) = dble(ix)/n1
              r(2) = dble(iy)/n2
              r(3) = dble(iz)/n3 ! in direct coordinate
              r(:) = matmul(box, r) ! convert to cartesian
              dr = r - RI
              ! PBC shifts
              ds = matmul(transpose(box_inv), dr)
              ds(1) = ds(1) - floor(ds(1)+0.5d0)
              ds(2) = ds(2) - floor(ds(2)+0.5d0)
              ds(3) = ds(3) - floor(ds(3)+0.5d0)

              ! three loops over neighbour cells for diffuse basis
              ! for most gaussian pair, these loops should do nothing
              do ixcell = ixlo, ixhi
              dsp(1) = ds(1) + ixcell
              do iycell = iylo, iyhi
              dsp(2) = ds(2) + iycell
              do izcell = izlo, izhi
              dsp(3) = ds(3) + izcell

              dr = matmul(box, dsp)
              foo = dot_product(dr,dr)*(alpha1+alpha2)
              if (foo > exp_thresh) then
                goto 111 ! skip point image i'
              else
                foo = volelm*dbl_mb(ivr+i-1)*fac*exp(-foo)
              endif
              xA = dr(1) + RI(1) - Apos(1)
              yA = dr(2) + RI(2) - Apos(2)
              zA = dr(3) + RI(3) - Apos(3)
              xRB = dr(1) + RI(1) - RBpos(1)
              yRB = dr(2) + RI(2) - RBpos(2)
              zRB = dr(3) + RI(3) - RBpos(3)
              ! Angular factors
              do a = 1, ncartA
                dbl_mb(iacart+a-1) = xA**int_mb(i1cart+a-1) 
     $                             * yA**int_mb(j1cart+a-1) 
     $                             * zA**int_mb(k1cart+a-1)
              end do
              do b = 1, ncartB
                dbl_mb(ibcart+b-1) = xRB**int_mb(i2cart+b-1)
     $                              *yRB**int_mb(j2cart+b-1)
     $                              *zRB**int_mb(k2cart+b-1)
              end do
              do b = 1, ncartB
                if (lsameshell .and. m==n) then
                  alo = b
                else
                  alo = 1
                endif
                do a = alo, ncartA
                  Final(m,n,a,b) = Final(m,n,a,b) 
     $                             + dbl_mb(iacart+a-1)*
     $                               dbl_mb(ibcart+b-1)*foo
                enddo
              enddo

111           continue
              enddo
              enddo
              enddo
             
333           continue 
              ! increment of real space indices
              ix = ix + 1
              if ( ix == n1 ) then
                ix = 0
                iy = iy + 1
                if ( iy == n2 ) then    
                  iy = 0
                  iz = iz + 1
                endif
              endif
            enddo ! end the loop over grid points
999         continue
          enddo ! end the loop over gaussian pairs
        enddo

        ! deal with the transpose symmetry
        !if ( lsameshell ) then
        !  do n = 1, nexpB
        !    do m = n, nexpA
        !      do b = 1, ncartB
        !        if ( m==n ) then
        !          alo = b+1
        !        else
        !          alo = 1
        !        endif
        !        do a = alo, ncartA
        !          Final(n,m,b,a) = Final(m,n,a,b)
        !        enddo
        !      enddo
        !    enddo
        !  enddo
        !endif

        if( .not. ma_chop_stack(li1cart) )
     $    call md_abort('uniint: Failed to deallocate scratch memory',0)

        return
      end subroutine 

      ! This sort out the order of the cartesians in a particular shell
      ! a stupid wrapper inherited from gamess_emb
      subroutine set_angular_index( l, icart, jcart, kcart, ncart )
      implicit none
        integer :: l, ncart, a, Nxyz(3)
        integer :: icart(ncart), jcart(ncart), kcart(ncart)

        do a = 1, ncart
          call getNxyz(l,a,Nxyz)
          icart(a) = Nxyz(1)
          jcart(a) = Nxyz(2)
          kcart(a) = Nxyz(3)
        enddo

        return
      end subroutine set_angular_index

      !*****************************************************
      ! set up the k-points information for the integration
      !*****************************************************
      subroutine setup_kpoints(box,box_inv, nx,ny,nz, npts, kx,ky,kz)
      implicit none
      include 'embcom.fh'
        double precision :: box(3,3), box_inv(3,3)
        integer :: nx, ny, nz, npts
        double precision :: kx(nx/2+1), ky(ny), kz(nz)
        integer :: ix, iy, iz, ipt, k1, k2, k3

        ipt = 1
        do iz = 1, nz
          k3 = iz - 1
          if ( k3 > nz/2 ) k3 = k3 - nz
          do iy = 1, ny
            k2 = iy - 1
            if ( k2 > ny/2 ) k2 = k2 - ny
            do ix = 1, nx/2+1
              k1 = ix - 1
              kx(ipt) = ( k1*box_inv(1,1) + k2*box_inv(1,2) +
     $                    k3*box_inv(1,3) )*2.0d0*PI
              ky(ipt) = ( k1*box_inv(2,1) + k2*box_inv(2,2) +
     $                    k3*box_inv(2,3) )*2.0d0*PI
              kz(ipt) = ( k1*box_inv(3,1) + k2*box_inv(3,2) +
     $                    k3*box_inv(3,3) )*2.0d0*PI
              ipt = ipt + 1
            enddo
          enddo
        enddo

        return
      end subroutine setup_kpoints

      !***********************************************
      ! calculate \int e^ikx * x^l * exp^(-alpha*x^2)
      !***********************************************
      subroutine int_pw_gauss(k, l, alpha, f_kla)
      implicit none
      include 'embcom.fh'
        double precision :: k, alpha
        integer :: l
        double complex :: f_kla(0:l)
        double complex, parameter :: iunit = cmplx(0,1)
        integer :: i

        f_kla(0) = cmplx( sqrt(PI/alpha)*exp(-k**2/4.d0/alpha), 0.0d0 )
        if ( l == 0 ) goto 666
        f_kla(1) = iunit*k/2/alpha * f_kla(0)
        if ( l == 1 ) goto 666

        do i = 2, l
          f_kla(i) = (i-1)/2.d0/alpha*f_kla(i-2) +
     $               iunit*k/2.d0/alpha*f_kla(i-1)
        enddo

666     continue

        return
      end subroutine int_pw_gauss

      !*************************************************
      ! setup the index bound for realspace integration
      !*************************************************
      subroutine set_ibound( rr, alpha, ilo, ihi, thresh )
      implicit none
      double precision :: rr, alpha, thresh
      integer :: ilo, ihi
        ilo = 0
        ihi = 0
        if ( rr*alpha < thresh ) then
          ilo = -1
          ihi = 1
        endif
        return
      end subroutine set_ibound


C============================================================
C This subroutine add the nonlocal Vemb matrix element to
C the one-electron integral matrix
C============================================================
      subroutine add_nlemb_elem(ish, i_basis, jsh, j_basis, i_gen, j_gen
     &                          , Li, Lj, Nints, V)
      implicit none
#include "embcom.fh"
#include "global.fh"
#include "stdio.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer ish,i_basis,jsh,j_basis,i_gen,j_gen,Nints,Li,Lj
      double precision V(Nints)
      integer ifirst,ilast,jfirst,jlast
      integer nbhandl, lbuff, ibuff, nbhandle, Li2, Lj2, i, j, nn, nnp

      if ( i_gen > 1 .or. j_gen > 1 ) then
        write(luout,*) 
     &       'NLVEMB WARNING: only segmented basis set is supported'
      endif
      if ( i_basis .ne. j_basis ) then
        write(luout,*) 
     &       'NLVEMB WARNING: multple basis sets were never tested'
      endif

      if ( .not. bas_cn2bfr(i_basis,ish,ifirst,ilast) ) 
     &   call errquit('add_nlemb_elem: fail to retrieve
     & basis function index', 0, BASIS_ERR )

      if ( .not. bas_cn2bfr(j_basis,jsh,jfirst,jlast) )
     &   call errquit('add_nlemb_elem: fail to retrieve
     & basis function index', 0, BASIS_ERR )

      ! allocate scratch buffer for potential matrix element
      if (.not.ma_push_get(mt_dbl,Nints,'nonloc Vemb scr',lbuff,ibuff)) 
     &   call errquit('add_nlemb_elem: fail to allocate buffer'
     &                ,0,MA_ERR)
      
      ! fetch nonlocal Vemb block, note take the transpose part,
      ! assuming the transpose symmetry for the V matrix. 
      ! the resulting buff array should be arranged by row, instead of
      ! column.
      call ga_nbget( g_pot, jfirst, jlast, ifirst, ilast, dbl_mb(ibuff)
     &              ,jlast-jfirst+1, nbhandle )
      ! modify the one-electron integrals
      Li2 = (Li+1)*(Li+2)/2
      Lj2 = (Lj+1)*(Lj+2)/2
      nn = 0
      do i = 1, Li2
        do j = 1, Lj2
          nn = nn + 1
          V(nn) = V(nn) + dbl_mb(ibuff+nn-1)
        enddo
      enddo

      ! deallocate memory
      if (.not. ma_pop_stack(lbuff)) 
     &   call errquit('add_nlemb_elem: fail to deallocate buffer'
     &                ,0,MA_ERR)

      return

      end subroutine add_nlemb_elem
