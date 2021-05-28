      module bertha_wrapper

      use, intrinsic :: iso_c_binding

      use spec, only: coord, znuc, cnuc, expset, iznuc, icrge, kvals,
     > nkappa, lmaxx, nfunct, ncent, ndim, nshift, nocc,
     > iter, iall, irun, readinput, verbose, dumpfiles
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

      implicit none

      integer :: intype, fitcase, testfit, igrid, iprint, mactve,
     +           maxit, nactve, ncore, nfitt, densitydiff
      real :: amass(ncentm), sfact1, sfact2, erep, thresh
      real :: eecoul, eexc 
      character(len=7) :: expo

      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_init(fnameinput, verbosity, dumpfilesin) 

      implicit none

      character*(*), intent(in) :: fnameinput

      integer (c_int) :: verbosity, dumpfilesin

      integer :: lammax

      verbose = verbosity
      dumpfiles = dumpfilesin

      if (verbose .eq. 1) call writeversion ()

      ndim = 0
      maxit = 2
      nshift = 0

      call factrl ()

      lammax = 4 * maxlqn + 10
      call gfinit(lammax)

      if(.not.readinput (0, intype, amass, ncore, mactve, nactve,
     +                   igrid, nopen, iprint, maxit, expo,
     +                   fitcase, testfit, sfact1, sfact2, fnameinput)) then
        call flush(6)
        return
      endif

      if (verbose .eq. 1) write(6,11) ndim
   11 format(2x,'total dimension of the matrix: ',i6/)
      
      call flush(6)

      end subroutine bertha_init

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_main(fittcoefffname, vctfilename, ovapfilename,
     +                       fittfname, eigen, ovap_ptr, eige_ptr, 
     +                       fock_ptr) 

      implicit none

      real (c_double) :: eigen(*)

      character*(*), intent(in) :: fittcoefffname, vctfilename,
     >   ovapfilename, fittfname
      real (c_double), intent(out) :: ovap_ptr(*), eige_ptr(*), 
     +                                fock_ptr(*)

      complex, parameter :: czero = (0.0,0.0)

      integer :: counter, icent, j, ii, jj, i, iter, 
     + nuca, nucb, nucc, nuccnt, rc, istate, a, iocc

      real, allocatable :: v(:), zz(:)

      real :: monoelectronenergy, centfun(ncentm), 
     + dist, r1x, r1y, r1z, d1, r2x, 
     + r2y, r2z, d2, scprod, ang, diff, eold, oldecoul, rthresh1, 
     + rthresh2, efittcoulplusexc, efittcoul, etemp, excfit, 
     + exc, cputime11, timetomt1, ffre,
     + ctimetomt1, timetomt2, ctimetomt2, wtimeiter, 
     + ctimeiter, ctime1, time1, ffim, ecoul

      real :: denswalltime, denscputime
      real :: energy1, wtimeio, ctimeio, energy2
      real :: coulenergy

      complex, allocatable :: olddensity(:,:)
      real :: maxdiff , percmaxdiff

      real, external :: walltime, cputime

      logical :: filefound

      cputime11 = cputime ()

      if (dumpfiles .eq. 1) then 
        open (unit=11, file=fittcoefffname, err=9999,
     +        status ="unknown")
      endif

      call getmem ()

      if (densitydiff .eq. 1) then 
        allocate(olddensity(ndim, ndim),stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        olddensity = (0.0, 0.0)
      endif

      call initijkvec()

      call initfitting(fitcase,nfitt,fittfname)
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

      if (verbose .eq. 1) then 
       write(6,*)
     &  '************************************************************'
       write(6,*)
       write(6,*)
       write(6,*)
     &  '        molecular geometry a: cartesian coordinates         '
       write(6,*)
       write(6,*)
     &  'center        x-coordinate    y-coordinate    z-coordinate'
       write(6,"(60('='))")
       do nuca = 1,ncent
         write(6,1118) elmnt(iznuc(nuca)),
     &                 coord(1,nuca),coord(2,nuca),coord(3,nuca)
 1118    format(3x,a,6x,f14.6,2x,f14.6,2x,f14.6)
         write(6,*)
     &   '------------------------------------------------------------'
       end do
       
       write(6,*)
       write(6,*)
       write(6,*)
     &  '                molecular geometry b: r-matrix              '
       write(6,*)
       write(6,*)
     &  'c1  c2       bond distance  c1  c2  c3     angle (degrees)  '
       write(6,"(60('='))")
      endif

      nuccnt = 1
      do nucb = 2,ncent
        r1x = coord(1,nucb)-coord(1,nuccnt)
        r1y = coord(2,nucb)-coord(2,nuccnt)
        r1z = coord(3,nucb)-coord(3,nuccnt)
        d1 = sqrt(r1x**2+r1y**2+r1z**2)
        if (verbose .eq. 1) write(6,1119) elmnt(iznuc(nuccnt)),elmnt(iznuc(nucb)),d1
        do nucc = 2,nucb-1
          r2x = coord(1,nucc)-coord(1,nuccnt)
          r2y = coord(2,nucc)-coord(2,nuccnt)
          r2z = coord(3,nucc)-coord(3,nuccnt)
          d2 = sqrt(r2x**2+r2y**2+r2z**2)
          scprod = (r1x*r2x+r1y*r2y+r1z*r2z)
          ang = acos(scprod/(d1*d2))*(180.0/pi)
          if (verbose .eq. 1) write(6,1120) elmnt(iznuc(nucb)),
     >                  elmnt(iznuc(nuccnt)),
     >                  elmnt(iznuc(nucc)),ang
        end do
        if (verbose .eq. 1) write(6,*)
     &  '------------------------------------------------------------'
      end do

 1119 format(2x,a,2x,a,4x,f14.6)
 1120 format(30x,a,2x,a,2x,a,2x,f14.6)

      if (verbose .eq. 1) then
       write(6,*)
       write(6,*)
       write(6,*)
       write(6,*)
      endif

      if (verbose .eq. -1) then 
        write(6,*) "Initialze the atomic densities on each center"
      endif

      if (irun.eq.1) then
        ioccm0 = nshift + 1
        if (verbose .eq. 1) then
         write(6,*)
     &    '          initialze the atomic densities on each center   '
         write(6,*)
         write(6,*)
        endif 

        do icent = 1,ncent
          call atomic(amass,iznuc,icrge,lmaxx,nkappa,icent,nshift)
        end do

        !open (unit=67, file="testc")
        !write(67, *) c
        !close(67)

        ioccm0 = ioccm0 - nshift
        call dens()
        if (densitydiff .eq. 1) then 
          if (allocated(olddensity)) olddensity = c
        endif


        ioccm0 = nocc
      else
        inquire(file=vctfilename, exist=filefound)
        if (.not.filefound) call errmsg('Restart file not found.')

        open (unit=10, file=vctfilename, err=9999, status='unknown')
        rewind(unit=10)
        read (10, '(I5,I5)') nocc, ndim
        do iocc = nshift+1,nshift+nocc
          read(10,'(E23.15)') eigen(iocc)
          do i=1,ndim
            read(10,'(E17.10,1X,E17.10)') c(i,iocc)
          end do
        end do
        close(unit=10)

        ioccm0 = nocc
        call dens()
        if (densitydiff .eq. 1) then 
          if (allocated(olddensity)) olddensity = c
        endif
      end if

      if (verbose .eq. -1) then 
        write(6,*) "Done"
      endif

      if (irun.eq.2) then
        ioccm0 = nocc
        allocate(v(nfitt),stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        allocate(zz(nfitt),stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        open(unit=11, file=fittcoefffname, status='unknown')
        rewind(11)
        read(11,*) v

        a=0
        excfit=0.0

        call dftinit_fit(excfit,igrid,expo,v,nfitt,fitcase,zz,a)

        v=v+zz

        fock=(0.0,0.0) 

        call jmatfitt(fock, ndim, v, nfitt, fitcase)

        deallocate(v)
        deallocate (zz)

        call oneel(monoelectronenergy)

        call diagonalizer (eigen)

        call dens()
        if (densitydiff .eq. 1) then 
          if (allocated(olddensity)) olddensity = c
        endif

      endif

      diff = 1.0
      eold = 1.0
      oldecoul = 1.0
      rthresh1 = 0.00001
      rthresh2 = 0.0000001

  777 format
     > ("=============================================================")
      do iter = 1,maxit
        if (verbose .eq. -1) then 
          write(6,*) "Iteration ", iter
        endif

        if (verbose .eq. 1) write (0, 777)
        wtimeiter = walltime ()
        ctimeiter = cputime ()
        fock = czero
        ovap = czero

        time1 = walltime ()
        ctime1 = cputime()
        call hden ()
        time1 = walltime () - time1
        ctime1 = cputime() - ctime1

        if (verbose .eq. 1) write (0, 6661 ) time1, ctime1
 6661   format ('Time for hden walltime                   ', f10.2,
     +  '  cputime ', f10.2,/)

        etotal = 0.0
        timetomt1 = walltime ()
        ctimetomt1 = cputime ()

        allocate(v(nfitt),stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        allocate(zz(nfitt), stat=rc)
        if (rc.ne.0) call errmsg('error allocating memory')
        call v_calculation(nfitt,fitcase,v)

        timetomt2 = walltime ()
        ctimetomt2 = cputime ()
        call factlu(nfitt,v)
        timetomt2 = walltime () - timetomt2
        ctimetomt2 = cputime () - ctimetomt2

        if (verbose .eq. 1) write (0, 5678) timetomt2, ctimetomt2
 5678   format ('Time for the factlu walltime             ', 
     +          f10.2, '  cputime ', f10.2)

        if (dumpfiles .eq. 1) write (11, *) v
        if (dumpfiles .eq. 1) rewind (11)

        call jfitenergy(v,nfitt,efittcoul)
        a=0
        excfit=0.0

        timetomt2 = walltime ()
        ctimetomt2 = cputime ()
        call dftinit_fit(excfit,igrid,expo,v,nfitt,fitcase,zz,a)
        timetomt2 = walltime () - timetomt2
        ctimetomt2 = cputime () - ctimetomt2

        if (verbose .eq. 1) write (0, 5680) timetomt2, ctimetomt2
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

        if (verbose .eq. 1) write (0, 5679) timetomt1, ctimetomt1
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

        if (verbose .eq. 1) then
         write(6,*) '******** energies values ********* '
         write(6,*) ' '
        endif

        if(testfit /= 0) then
          if (verbose .eq. 1) 
     >      write(6,*) '******** ATTENTION rjmat is on. This should be used only for benchmark********* '
          call rjmat(ndim,c,coulenergy) 
          if (verbose .eq. 1) then
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
          endif
        else
          if (verbose .eq. 1) then
           write(6,*) '******** energies fitted values ********* '
           write(6,*) ' '
           write(6,*) 'coulomb  energy(fitting)     : ',efittcoul/2.0d0
           write(6,*) 'Exc(rhofit)                  : ',excfit
           write(6,*) ' '
          endif
        endif

        eecoul  = efittcoul
        eexc    = excfit


        energy2 = excfit - etemp

        time1 = walltime ()
        ctime1 = cputime ()
        call oneel(monoelectronenergy)
        time1 = walltime () - time1
        ctime1 = cputime () - ctime1

        if (dumpfiles .eq. 1) then
         open (unit=99, file=ovapfilename, err=9999, status='unknown')
         write (99, '(I5)') ndim
         do i=1,ndim
           do j=1,ndim
             write(99,'(E22.10E4,1X,E22.10E4)') ovap(i, j)
           end do
         end do
         close(unit=99)
        endif

        counter = 1
        do i = 1, ndim
          do j = 1, ndim
c            write (6, '("(",F10.5," ",F10.5,"i)")', advance="no") ovap(i, j)
c            write (6, '("(",F10.5," ",F10.5,"i)")') ovap(i, j)
             ovap_ptr(counter) = real(ovap(i,j))
             ovap_ptr(counter+1) = imag(ovap(i,j))
             counter = counter + 2
          enddo
c          write (6, *) " "
        enddo
        
c        ovap_ptr = c_loc(ovap(1,1))

        if (verbose .eq. 1) write (0, 6664 ) time1, ctime1
 6664   format ('Time for oneel walltime                  ', f10.2,
     +          '  cputime ', f10.2, /)
        time1 = walltime ()
        ctime1 = cputime ()


c        if (((iter.gt.1).and.(sfact.ne.(0.0))).or.(irun.eq.0).
        if (((iter.gt.1).and.(abs(sfact).gt.(1.0d-17))).or.(irun.eq.0).
     >      or.(irun.eq.2))
     >     call shftlv()
        time1 = walltime () - time1
        ctime1 = cputime () - ctime1
        if (verbose .eq. 1) write (0, 5564 ) time1, ctime1
 5564   format ('Time for shftlv                          ', f10.2,
     +          '  cputime ', f10.2, /)

        do i = 1,ndim
          ffre = real(fock(i,i))
          fock(i,i) = cmplx(ffre,0.0)
        end do

        counter = 1
        do i = 1, ndim
          do j = 1, ndim
             fock_ptr(counter) = real(fock(i,j))
             fock_ptr(counter+1) = imag(fock(i,j))
             counter = counter + 2
          enddo
        enddo

#ifdef DUMPFOCKMTX
        open(unit=44, file="fockmtx.txt")
        do i = 1,ndim
          do j = 1,ndim
            write (44, '(I5,1X,I5,E22.10E4,1X,E22.10E4)' ) i, j, fock(i,j)
          end do
        end do
        close(44)
#endif

        call diagonalizer (eigen)

        wtimeio = walltime ()
        ctimeio = cputime ()

        if (dumpfiles .eq. 1) then
         open (unit=10, file=vctfilename, err=9999, status='unknown')
         rewind(unit=10)

         write (10, '(I5,I5)') nocc, ndim
         do iocc = nshift+1,nshift+nocc
           write(10,'(E23.15)') eigen(iocc)
           do i=1,ndim
             write(10,'(E17.10,1X,E17.10)') c(i,iocc)
           end do
         end do

c        write (10, '(I5,I5)') ndim, ndim
c        do iocc = 1,ndim
c          write(10,'(E18.10)') eigen(iocc)
c          do i=1,ndim
c            write(10,'(E17.10,1X,E17.10)') c(i,iocc)
c          end do
c        end do
 
         close(unit=10)
        endif

        counter = 1
        do j = 1, ndim
          do i = 1, ndim
             eige_ptr(counter) = real(c(j,i))
             eige_ptr(counter+1) = imag(c(j,i))
             counter = counter + 2
          enddo
        enddo
        
        energy1 = 0.0
        do istate = 1,nocc
          energy1 = energy1 + eigen(istate+nshift)
        end do
        etotal = energy1 + energy2
        
        wtimeio = walltime () - wtimeio
        ctimeio = cputime () - ctimeio
        
        if (verbose .eq. 1) write (0, 3567) wtimeio, ctimeio
 3567   format ("Time for the io walltime                 ", 
     +        f10.2, "  cputime ", f10.2, /)

        if (iprint.eq.2) then
          if (verbose .eq. 1) then 
           write(6,*) 'iteration number: ',iter
           write(6,*) 'level shift parameter = ',sfact
           do istate = 1,nocc+nopen
             write(6,*) 'eigenvalue number, e:',
     &                  istate,eigen(istate+nshift)-sfact
           end do
           write(6,*) 'lumo               e:',
     &                istate,eigen(nocc+nopen+nshift+1)
           write(6,*) 'total electronic energy = ',etotal-(sfact*nocc)
           write(6,*) 'nuclear repulsion energy = ',erep
           write(6,*) 'total energy = ',etotal+erep-(sfact*nocc)
           write(6,*) ' '
           write(6,*) '-----------------------------------------------'
          endif
        end if

        if (verbose .eq. -1) then
           write(6,'(A, " ", F20.10)') '  total electronic energy   = ',etotal-(sfact*nocc)
           write(6,'(A, " ", F20.10)') '  nuclear repulsion energy  = ',erep
           write(6,'(A, " ", F20.10)') '  total energy              = ',etotal+erep-(sfact*nocc)
        endif

        diff = abs(eold-etotal)/(abs(etotal)+1.0)
        if (diff.lt.thresh.and.iall.eq.2) then
          goto 2000
        else
          continue
        endif

        denswalltime = walltime ()
        denscputime = cputime ()
        
        call dens()
        
        if (densitydiff .eq. 1) then 
          maxdiff = 0.0 
          percmaxdiff = 0.0

          do j=1,ndim
            do i=1,ndim
              if (abs(olddensity(i, j) - c(i, j)) .gt. maxdiff) then 
                maxdiff = abs(olddensity(i, j) - c(i, j))
                percmaxdiff = 100.0 * (maxdiff/(abs(olddensity(i, j) + 
     +              c(i,j))/2.0))
              endif
            enddo
          enddo 

          if (verbose .eq. -1) then
             write(6,'(A, " ", E20.10, " (", F8.3, "%) ")') 
     +       '  maxdiff in density matrix = ',maxdiff, percmaxdiff
          endif

          if (allocated(olddensity)) olddensity = c
        endif

        denswalltime = walltime () - denswalltime
        denscputime = cputime () - denscputime

        if (verbose .eq. 1) write (0, 3457) denswalltime, denscputime
 3457   format ('Time for dens walltime                   ', f10.2,
     +          '  cputime ', f10.2, /)

        if (diff.lt.rthresh1.and.iall.eq.0) then
          iall = 1
          sfact = sfact1
        else if (diff.lt.rthresh2.and.iall.eq.1) then
          iall = 2
          sfact = sfact2
        end if
        eold = etotal

        wtimeiter = walltime () - wtimeiter
        ctimeiter = cputime () - ctimeiter

        if (verbose .eq. 1) write (0, 777)
        if (verbose .eq. 1) write (0, 6666) iter, wtimeiter, ctimeiter
 6666   format ("Total time for iteration ", I4, " walltime   ",
     +          F10.2, "  cputime ", F10.2, /)
      end do

      if (verbose .eq. 1) write(6,*) 'convergence not obtained in the maximum number'
      if (verbose .eq. 1) write(6,*) 'of iterations'

 2000 continue
      if (iter.le.maxit) then
        wtimeiter = walltime () - wtimeiter
        ctimeiter = cputime () - ctimeiter

        if (verbose .eq. 1) write (0, 777)
        if (verbose .eq. 1) write (0, 6666) iter, wtimeiter, ctimeiter
      endif

      cputime11 = cputime11 - cputime ()

      if (verbose .eq. 1) then
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
      endif

      iall = 2
      irun = 0

c     call getresourcesusage ()  

      if (densitydiff .eq. 1) then 
        deallocate(olddensity)
      endif

      wtimeiter = walltime ()
      !call gencubefile(ncent,coord,iznuc,nfitt,fitcase)
      wtimeiter = walltime () - wtimeiter

      if (verbose .eq. 1) write (0, 6679) wtimeiter
 6679 format ('Time for cube is: ', f8.4)

      call flush(6)

      return 

 9999 call errmsg("Error open file")

      end subroutine 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_finalize()

      call fremem()

      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_realtime_init() 

      use realtime

      implicit none

      call allocate_real_time(ndim)

      end subroutine bertha_realtime_init

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_realtime_finalize() 

      use realtime

      implicit none

      call deallocate_real_time()

      end subroutine bertha_realtime_finalize

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_realtime_dipolematrix (direction, norm, 
     +                                         vext_ptr)

      use realtime

      implicit none

      integer, intent(in) :: direction, norm
      real (c_double), intent(out) :: vext_ptr(*) 

      integer :: counter, i, j
      
      logical :: normalise

      if (norm.eq.1) normalise = .true.
      if (norm.eq.0) normalise = .false.

      ! come se fosse px a regime avremo px py pz e quindi
      ! tre matrici  vext (2 = z, 3 = y, 4 = x)
      call dipolematrix (direction, normalise)

      counter = 1
      do i = 1, ndim
        do j = 1, ndim
           vext_ptr(counter) = real(vext(i,j))
           vext_ptr(counter+1) = imag(vext(i,j))
           counter = counter + 2
        enddo
      enddo

      end subroutine bertha_realtime_dipolematrix 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_eps (x, y, z, val)

      implicit none

      real (c_double), intent(in) :: x, y, z
      real (c_double), intent(out) :: val

      call eps (x, y, z, val)

      end subroutine bertha_eps

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_realtime_fock (dens_ptr, fock_ptr)

      use fmat, only: fock
      use coeff, only: c
      use realtime

      implicit none 

      real (c_double), intent(out) :: fock_ptr(*)
      real (c_double), intent(in) :: dens_ptr(*) 

      real :: excfit, efittcoul, monoelectronenergy, ffre
      real :: efittcoulplusexc, traceDFock
      complex, parameter :: czero = (0.0,0.0)
      real, allocatable :: v(:), zz(:)
      integer :: counter, i, j, a, rc

      counter = 1
      do i = 1, ndim
        do j = 1, ndim
           c(i, j) = cmplx(dens_ptr(counter), dens_ptr(counter+1))
           counter = counter + 2
        enddo
      enddo

      fock = czero
      call hden ()
      allocate(v(nfitt),stat=rc)
      if (rc.ne.0) call errmsg('error allocating memory')
      allocate(zz(nfitt), stat=rc)
      if (rc.ne.0) call errmsg('error allocating memory')
      call v_calculation(nfitt,fitcase,v)
      call factlu(nfitt,v)
      call jfitenergy(v,nfitt,efittcoul)
      a = 0
      excfit = 0.0
      call dftinit_fit(excfit,igrid,expo,v,nfitt,fitcase,zz,a)
      v = zz + v
      call jmatfitt(fock,ndim,v,nfitt,fitcase)
      deallocate(v)
      deallocate(zz)

      do i = 1,nshift
        do j = 1,i
          fock(j,i) = conjg(fock(i,j))
          fock(j+nshift,i+nshift) = conjg(fock(i+nshift,j+nshift))
        end do
      end do

! #### Calculation of Tr(D(V_xc+J))    
 
      efittcoulplusexc = 0.0
      do i = 1,ndim
          do j = 1,ndim
              efittcoulplusexc = efittcoulplusexc +
     +        real(c(i,j)*fock(i,j))
          end do
      end do

! #### end  Calculation of Tr(D(V_xc+J))    

      call oneel(monoelectronenergy)

      do i = 1,nshift
        do j = 1,i
          fock(j,i) = conjg(fock(i,j))
          fock(j+nshift,i+nshift) = conjg(fock(i+nshift,j+nshift))
        end do
      end do

      do i = 1,ndim
        ffre = real(fock(i,i))
        fock(i,i) = cmplx(ffre,0.0)
      end do

! #### Calculation of Tr(D Fock)
      traceDFock = 0.0
      do i = 1,ndim
          do j = 1,ndim
              traceDFock =  traceDFock +
     +        real(c(i,j)*fock(i,j))
          end do
      end do
! #### end Calculation of Tr(D Fock)
! 
! ####  Calculation of total energy given a density matrix

      etotal = traceDFock - efittcoulplusexc + 0.5d0*efittcoul + excfit 
      eecoul  = efittcoul 
      eexc    = excfit 

! ####  end Calculation of total energy given a density matrix

      counter = 1
      do i = 1, ndim
        do j = 1, ndim
           fock_ptr(counter) = real(fock(i,j))
           fock_ptr(counter+1) = imag(fock(i,j))
           counter = counter + 2
        enddo
      enddo

      end subroutine bertha_realtime_fock

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_density_to_cube (dens_ptr, margin, drx, dry, 
     >    drz, fname, fittfname)

      implicit none

      character*(*), intent(in) :: fname, fittfname
      real (c_double), intent(in) :: dens_ptr(*) 
      real, intent(in) :: drx, dry, drz, margin

      real :: centfun(ncentm)
      integer :: counter, i, j
      real, allocatable :: v(:), zz(:)

      real, dimension (3) :: ri, rf, dr
      integer, dimension (3) :: nr

      counter = 1
      do i = 1, ndim
        do j = 1, ndim
           c(i, j) = cmplx(dens_ptr(counter), dens_ptr(counter+1))
           counter = counter + 2
        enddo
      enddo

      ri(1) = coord(1, 1)
      ri(2) = coord(2, 1)
      ri(3) = coord(3, 1)

      rf(1) = coord(1, 1)
      rf(2) = coord(2, 1)
      rf(3) = coord(3, 1)

      do i=1,ncent
        !print *, coord(1, i), coord(2, i), coord(3, i)

        if (coord(1, i) .lt. ri(1)) ri(1) = coord(1, i) 
        if (coord(2, i) .lt. ri(2)) ri(2) = coord(2, i)
        if (coord(3, i) .lt. ri(3)) ri(3) = coord(3, i)

        if (coord(1, i) .gt. rf(1)) rf(1) = coord(1, i) 
        if (coord(2, i) .gt. rf(2)) rf(2) = coord(2, i)
        if (coord(3, i) .gt. rf(3)) rf(3) = coord(3, i)
 
      enddo

      !print *, drx, dry, drz
      !print *, "Min: ", ri(1), ri(2), ri(3)
      !print *, "Max: ", rf(1), rf(2), rf(3)

      dr(1) = drx
      dr(2) = dry 
      dr(3) = drz

      !print *, "dr: ", dr(1), dr(2), dr(3)

      ri(1) = ri(1) - margin
      ri(2) = ri(2) - margin
      ri(3) = ri(3) - margin

      rf(1) = rf(1) + margin
      rf(2) = rf(2) + margin
      rf(3) = rf(3) + margin

      nr = INT( (rf - ri) / dr + 1.0 )

      call initijkvec()
      call initfitting(fitcase, nfitt, fittfname)
      call flabel(centfun)
      allocate (v(nfitt),zz(nfitt))

      ! rho
      call hden()
      call v_calculation(nfitt,fitcase,v)
      call factlu(nfitt,v)
      ! write(filename, '(A)') "fit-rho.cube"
      ! print *, "Writing cubefile: ", filename
      call cubefile (ncent,coord,iznuc,v,nfitt,fitcase,
     +               ri,nr,dr,fname) 

      !deallocate(denh)
      deallocate(zz)
      deallocate(v) 

      end subroutine bertha_density_to_cube 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_get_coord (i, x, y, z, an)

      implicit none

      real, intent(out) :: x, y, z, an
      integer, intent(in) :: i

      !print *, "Here"

      if ( i .le. ncent) then
        x = coord(1, i)
        y = coord(2, i)
        z = coord(3, i)

        an = znuc(i)
      endif

      end subroutine bertha_get_coord

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bertha_density_to_cube_limit (dens_ptr, ri1, ri2, ri3, 
     >    rf1, rf2, rf3, drx, dry, drz, fname, fittfname)

      implicit none

      character*(*), intent(in) :: fname, fittfname
      real (c_double), intent(in) :: dens_ptr(*) 
      real, intent(in) :: drx, dry, drz, ri1, ri2, ri3, rf1, rf2, rf3

      real :: centfun(ncentm)
      integer :: counter, i, j
      real, allocatable :: v(:), zz(:)

      real, dimension (3) :: ri, rf, dr
      integer, dimension (3) :: nr

      counter = 1
      do i = 1, ndim
        do j = 1, ndim
           c(i, j) = cmplx(dens_ptr(counter), dens_ptr(counter+1))
           counter = counter + 2
        enddo
      enddo

      ri(1) = ri1
      ri(2) = ri2
      ri(3) = ri3

      rf(1) = rf1
      rf(2) = rf2
      rf(3) = rf3

      dr(1) = drx
      dr(2) = dry 
      dr(3) = drz

      nr = INT( (rf - ri) / dr + 1.0 )

      call initijkvec()
      call initfitting(fitcase, nfitt, fittfname)
      call flabel(centfun)
      allocate (v(nfitt),zz(nfitt))

      call hden()
      call v_calculation(nfitt,fitcase,v)
      call factlu(nfitt,v)
      call cubefile (ncent,coord,iznuc,v,nfitt,fitcase,
     +               ri,nr,dr,fname) 

      deallocate(zz)
      deallocate(v) 

      end subroutine bertha_density_to_cube_limit

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end module
