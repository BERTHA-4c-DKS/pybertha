      module bertha_wrapper

      implicit none

      contains

      subroutine bertha_main()

      use spec, only: coord, znuc, cnuc, expset, iznuc, icrge, kvals,
     > nkappa, lmaxx, nfunct, ncent, ndim, nshift, nocc,
     > iter, iall, irun, readinput
      use fmat, only: fock, ovap
      use coeff, only: c
      use dcoeff, only: d
      use clabel, only: ioccm0
      use limits, only: ncentm, maxmv2, maxlqn, maxb, nkapm
      use constants, only: pi, bohr, elmnt
      use energy, only: etotal
      use shiftr, only: sfact
      use opensh, only: nopen
      use indsys
!     use timetot, only: tot
!     use looptime, only: loop3

      implicit none

      interface
        subroutine diagonalizer (eigen)
          real, allocatable, intent(inout) :: eigen(:)
        end subroutine diagonalizer
      end interface

      interface
        subroutine diagonalizer_iter (eigen)
          real, allocatable, intent(inout) :: eigen(:)
        end subroutine diagonalizer_iter
      end interface

      real, parameter :: tresh = 1.0e-11
      complex, parameter :: czero = (0.0,0.0)

      integer :: nfitt, fitcase, testfit, lammax, intype,
     + icent, j, ii, jj, i, iprint, ncore, mactve, nactve, 
     + igrid, nuca, nucb, nucc, nuccnt, rc, istate, a, iocc

      real, allocatable :: eigen(:), v(:), zz(:)

      real :: amass(ncentm), monoelectronenergy, centfun(ncentm), 
     + sfact1, sfact2, erep, dist, r1x, r1y, r1z, d1, r2x, 
     + r2y, r2z, d2, scprod, ang, diff, eold, oldecoul, rtresh1, 
     + rtresh2, efittcoulplusexc, efittcoul, etemp, excfit, 
     + exc, cputime11, timetomt1, ffre,
     + ctimetomt1, timetomt2, ctimetomt2, wtimeiter, 
     + ctimeiter, ctime1, time1, ffim, ecoul

      real :: denswalltime, denscputime
      real :: energy1, wtimeio, ctimeio, energy2
      real :: coulenergy

      character(len=7) :: expo
      character(len=11) :: vctfilename
      character(len=11) :: ovapfilename

      real, external :: walltime, cputime

      integer :: maxit

      logical :: filefound

      maxit = 2

      cputime11 = cputime ()

      call writeversion ()

      ! open (unit=98, file = "out.val", ERR=9999, status = "unknown")

      ndim = 0
      nshift = 0

      open (unit=11, file="fittcoeff.restart", err=9999,
     +      status ="unknown")

      write(vctfilename,'(a7)') "vct.out"
      write(ovapfilename,'(a8)') "ovap.out"
  
      call factrl ()

      lammax = 4 * maxlqn + 10
      call gfinit(lammax)

      if(.not.readinput (0, intype, amass, ncore, mactve, nactve,
     +                   igrid, nopen, iprint, maxit, expo,
     +                   fitcase, testfit, sfact1, sfact2)) then
        stop
      endif
      
      write(6,11) ndim
   11 format(2x,'total dimension of the matrix: ',i6/)

      call getmem ()
      allocate(eigen(ndim), stat=rc)

      call initijkvec()
      call initfitting(fitcase,nfitt)
      call flabel(centfun)
! nuclear repulsion energy
      erep = 0.0
      do nuca = 1,ncent
        do nucb = 1,nuca-1
          dist = sqrt((coord(1,nuca)-coord(1,nucb))**2 +
     + (coord(2,nuca)-coord(2,nucb))**2 +
     + (coord(3,nuca)-coord(3,nucb))**2)
          erep = erep + znuc(nuca) * znuc(nucb) / dist
        end do
      end do

      write(6,*)
     & '************************************************************'
      write(6,*)
      write(6,*)
      write(6,*)
     & '        molecular geometry a: cartesian coordinates         '
      write(6,*)
      write(6,*)
     & 'centre        x-coordinate    y-coordinate    z-coordinate'
      write(6,"(60('='))")
      do nuca = 1,ncent
        write(6,1118) elmnt(iznuc(nuca)),
     &                coord(1,nuca),coord(2,nuca),coord(3,nuca)
 1118   format(3x,a,6x,f14.6,2x,f14.6,2x,f14.6)
        write(6,*)
     &  '------------------------------------------------------------'
      end do

      write(6,*)
      write(6,*)
      write(6,*)
     & '                molecular geometry b: r-matrix              '
      write(6,*)
      write(6,*)
     & 'c1  c2       bond distance  c1  c2  c3     angle (degrees)  '
      write(6,"(60('='))")

      nuccnt = 1
      do nucb = 2,ncent
        r1x = coord(1,nucb)-coord(1,nuccnt)
        r1y = coord(2,nucb)-coord(2,nuccnt)
        r1z = coord(3,nucb)-coord(3,nuccnt)
        d1 = sqrt(r1x**2+r1y**2+r1z**2)
        write(6,1119) elmnt(iznuc(nuccnt)),elmnt(iznuc(nucb)),d1
        do nucc = 2,nucb-1
          r2x = coord(1,nucc)-coord(1,nuccnt)
          r2y = coord(2,nucc)-coord(2,nuccnt)
          r2z = coord(3,nucc)-coord(3,nuccnt)
          d2 = sqrt(r2x**2+r2y**2+r2z**2)
          scprod = (r1x*r2x+r1y*r2y+r1z*r2z)
          ang = acos(scprod/(d1*d2))*(180.0/pi)
          write(6,1120) elmnt(iznuc(nucb)),
     >                  elmnt(iznuc(nuccnt)),
     >                  elmnt(iznuc(nucc)),ang
        end do
        write(6,*)
     &  '------------------------------------------------------------'
      end do

 1119 format(2x,a,2x,a,4x,f14.6)
 1120 format(30x,a,2x,a,2x,a,2x,f14.6)

      write(6,*)
      write(6,*)
      write(6,*)
      write(6,*)

      if (irun.eq.1) then
        ioccm0 = nshift + 1
        write(6,*)
     &   '          initialze the atomic densities on each centre   '
        write(6,*)
        write(6,*)
        do icent = 1,ncent
          call atomic(amass,iznuc,icrge,lmaxx,nkappa,icent,nshift)
        end do

        !open (unit=67, file="testc")
        !write(67, *) c
        !close(67)

        ioccm0 = ioccm0 - nshift
        call dens()
        ioccm0 = nocc
      else
        inquire(file=vctfilename, exist=filefound)
        if (.not.filefound) call errmsg('Restart file not found.')

        open (unit=10, file=vctfilename, err=9999, status='unknown')
        rewind(unit=10)
        read (10, '(I5,I5)') nocc, ndim
        do iocc = nshift+1,nshift+nocc
          read(10,'(E18.10)') eigen(iocc)
          do i=1,ndim
            read(10,'(E17.10,1X,E17.10)') c(i,iocc)
          end do
        end do
        close(unit=10)

        ioccm0 = nocc
        call dens()
      end if

      if (irun.eq.2) then
        ioccm0 = nocc
        allocate(v(nfitt),zz(nfitt),stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        open(unit=11, file="fittcoeff.restart", status='unknown')
        rewind(11)
        read(11,*) v

        a=0
        excfit=0.0

        call dftinit_fit(excfit,igrid,expo,v,nfitt,fitcase,zz,a,.true.)

        v=v+zz

        fock=(0.0,0.0) 

        call jmatfitt(fock, ndim, v, nfitt, fitcase)

        deallocate(v)
        deallocate (zz)

        call oneel(monoelectronenergy)

        call diagonalizer (eigen)

        call dens()

      endif

      diff = 1.0
      eold = 1.0
      oldecoul = 1.0
      rtresh1 = 0.00001
      rtresh2 = 0.0000001

  777 format
     > ("=============================================================")
      do iter = 1,maxit
        write (0, 777)
        wtimeiter = walltime ()
        ctimeiter = cputime ()
        fock = czero
        ovap = czero

        time1 = walltime ()
        ctime1 = cputime()
        call hden ()
        time1 = walltime () - time1
        ctime1 = cputime() - ctime1

        write (0, 6661 ) time1, ctime1
 6661   format ('Time for hden walltime                   ', f10.2,
     +  '  cputime ', f10.2,/)

        etotal = 0.0
        timetomt1 = walltime ()
        ctimetomt1 = cputime ()

        allocate(v(nfitt),zz(nfitt), stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        call v_calculation(nfitt,fitcase,v)

        timetomt2 = walltime ()
        ctimetomt2 = cputime ()
        call factlu(nfitt,v)
        timetomt2 = walltime () - timetomt2
        ctimetomt2 = cputime () - ctimetomt2

        write (0, 5678) timetomt2, ctimetomt2
 5678   format ('Time for the factlu walltime             ', 
     +          f10.2, '  cputime ', f10.2)

        !write (11, *) v
        rewind (11)

        call jfitenergy(v,nfitt,efittcoul)
        a=0
        excfit=0.0

        timetomt2 = walltime ()
        ctimetomt2 = cputime ()
        call dftinit_fit(excfit,igrid,expo,v,nfitt,fitcase,zz,a,.true.)
        timetomt2 = walltime () - timetomt2
        ctimetomt2 = cputime () - ctimetomt2

        write (0, 5680) timetomt2, ctimetomt2
 5680   format ('Time for the zz walltime                 ', f10.2, 
     +          '  cputime ', f10.2)

        v = zz + v
        call jmatfitt(fock,ndim,v,nfitt,fitcase)
        deallocate(v)
        deallocate(zz)
        efittcoulplusexc = 0.0

        do ii = 1,ndim
          do jj = 1,ndim
            efittcoulplusexc = efittcoulplusexc + 
     +        real(c(ii,jj)*fock(ii,jj))
          end do
        end do
        etemp = efittcoulplusexc - efittcoul / 2
        timetomt1 = walltime () - timetomt1
        ctimetomt1 = cputime () - ctimetomt1

        write (0, 5679) timetomt1, ctimetomt1
 5679   format ('Total time for the fitting walltime      ', f10.2, 
     +          '  cputime ', f10.2, / )

        do i = 1,nshift
          do j = 1,i
            fock(j,i) = conjg(fock(i,j))
            fock(j+nshift,i+nshift) = conjg(fock(i+nshift,j+nshift))
          end do
        end do

        ecoul = 0.0
        do i = 1,ndim
          do j = 1,ndim
            ffre = real(fock(i,j))
            ffim = imag(fock(i,j))
            if (abs(ffre).lt.0.0000000001) ffre = 0.0
            if (abs(ffim).lt.0.0000000001) ffim = 0.0
            fock(i,j) = cmplx(ffre,ffim)
            ecoul = ecoul + real(c(i,j)*fock(i,j))
          end do
        end do
        ecoul = ecoul / 2
        exc = 0.0

        write(6,*) '******** energies values ********* '
        write(6,*) ' '

        if(testfit /= 0) then
          write(6,*) '******** ATTENTION rjmat is on. This should be used only for benchmark********* '
          call rjmat(ndim,c,coulenergy) 
          write(6,*) 'coulomb  energy              : ',coulenergy/2.0d0
          write(6,*) ' '
          write(6,*) '******** energies fitted values ********* '
          write(6,*) ' '
          write(6,*) 'coulomb  energy(fitting)     : ',efittcoul/2.0d0
          write(6,*) 'Exc(rhofit)                  : ',excfit
          write(6,*) ' '
          write(6,*) 'Error in Coulomb energy due to fitting: ',efittcoul/2.0d0-coulenergy/2.0
          write(6,*) ' '
          write(6,*) '******************************************************************************* '
          write(6,*) ' '
        else
          write(6,*) '******** energies fitted values ********* '
          write(6,*) ' '
          write(6,*) 'coulomb  energy(fitting)     : ',efittcoul/2.0d0
          write(6,*) 'Exc(rhofit)                  : ',excfit
          write(6,*) ' '
        endif

        energy2 = excfit - etemp

        time1 = walltime ()
        ctime1 = cputime ()
        call oneel(monoelectronenergy)
        time1 = walltime () - time1
        ctime1 = cputime () - ctime1

        open (unit=99, file=ovapfilename, err=9999, status='unknown')
        write (99, '(I5)') ndim
        do i=1,ndim
          do j=1,ndim
            write(99,'(E22.10E4,1X,E22.10E4)') ovap(i, j)
          end do
        end do
        close(unit=99)

        write (0, 6664 ) time1, ctime1
 6664   format ('Time for oneel walltime                  ', f10.2,
     +          '  cputime ', f10.2, /)
        time1 = walltime ()
        ctime1 = cputime ()

        if ((iter.gt.1.and.sfact.ne.0.0).or.irun.eq.0.or.irun.eq.2)
     >     call shftlv()
        time1 = walltime () - time1
        ctime1 = cputime () - ctime1
        write (0, 5564 ) time1, ctime1
 5564   format ('Time for shftlv                          ', f10.2,
     +          '  cputime ', f10.2, /)

        do i = 1,ndim
          ffre = real(fock(i,i))
          fock(i,i) = cmplx(ffre,0.0)
        end do

        call diagonalizer (eigen)

        wtimeio = walltime ()
        ctimeio = cputime ()

        open (unit=10, file=vctfilename, err=9999, status='unknown')
        rewind(unit=10)
        write (10, '(I5,I5)') nocc, ndim
        do iocc = nshift+1,nshift+nocc
          write(10,'(E18.10)') eigen(iocc)
          do i=1,ndim
            write(10,'(E17.10,1X,E17.10)') c(i,iocc)
          end do
        end do
        close(unit=10)
        
        energy1 = 0.0
        do istate = 1,nocc
          energy1 = energy1 + eigen(istate+nshift)
        end do
        etotal = energy1 + energy2
        
        wtimeio = walltime () - wtimeio
        ctimeio = cputime () - ctimeio
        
        write (0, 3567) wtimeio, ctimeio
 3567   format ("Time for the io walltime                 ", 
     +        f10.2, "  cputime ", f10.2, /)

        if (iprint.eq.2) then
          write(6,*) 'iteration number: ',iter
          write(6,*) 'level shift parameter = ',sfact
          do istate = 1,nocc+nopen
            write(6,*) 'eigenvalue number, e:',
     &                 istate,eigen(istate+nshift)-sfact
          end do
          write(6,*) 'lumo               e:',
     &               istate,eigen(nocc+nopen+nshift+1)
          write(6,*) 'total electronic energy = ',etotal-(sfact*nocc)
          write(6,*) 'nuclear repulsion energy = ',erep
          write(6,*) 'total energy = ',etotal+erep-(sfact*nocc)
          write(6,*) ' '
          write(6,*) '-----------------------------------------------'
        end if

        diff = abs(eold-etotal)/(abs(etotal)+1.0)
        if (diff.lt.tresh.and.iall.eq.2) then
          goto 2000
        else
          continue
        endif

        denswalltime = walltime ()
        denscputime = cputime ()
        call dens()
        denswalltime = walltime () - denswalltime
        denscputime = cputime () - denscputime

        write (0, 3457) denswalltime, denscputime
 3457   format ('Time for dens walltime                   ', f10.2,
     +          '  cputime ', f10.2, /)

        if (diff.lt.rtresh1.and.iall.eq.0) then
          iall = 1
          sfact = sfact1
        else if (diff.lt.rtresh2.and.iall.eq.1) then
          iall = 2
          sfact = sfact2
        end if
        eold = etotal

        wtimeiter = walltime () - wtimeiter
        ctimeiter = cputime () - ctimeiter

        write (0, 777)
        write (0, 6666) iter, wtimeiter, ctimeiter
 6666   format ("Total time for iteration ", I4, " walltime   ",
     +          F10.2, "  cputime ", F10.2, /)
      end do

      write(6,*) 'convergence not obtained in the maximum number'
      write(6,*) 'of iterations'

 2000 continue
      if (iter.le.maxit) then
        wtimeiter = walltime () - wtimeiter
        ctimeiter = cputime () - ctimeiter

        write (0, 777)
        write (0, 6666) iter, wtimeiter, ctimeiter
      endif

      cputime11 = cputime11 - cputime ()

      write(6,*) ' '
      write(6,*) '**********************************************'
      write(6,*) 'final output:'
      write(6,*) 'convergence obtained in',iter,'     iterations'
      write(6,*) ' '
      write(6,*) 'bond length = ',coord(3,2)
      write(6,*) 'total energy = ',etotal+erep-(nocc*sfact)
      write(6,*) ' '
      write(6,*) 'total cpu-time used (sec):', cputime11
      write(6,*) '**********************************************'

      iall = 2
      irun = 0

c     call getresourcesusage ()

      wtimeiter = walltime ()
      !call gencubefile(ncent,coord,iznuc,nfitt,fitcase)
      wtimeiter = walltime () - wtimeiter

      write (0, 6679) wtimeiter
 6679 format ('Time for cube is: ', f8.4)

      deallocate(eigen)
      call fremem()

      stop

 9999 call errmsg("Error open file")

      end subroutine 

      end module
