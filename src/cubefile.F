      subroutine cubefile (nnuc,coord,iznuc,v,nfitt,fitcase,
     +                     ri, nr, dr, filename)

      implicit none

      character(len=20),intent(in)                  :: filename
      integer,intent(in)                            :: nnuc
      double precision,dimension(3,nnuc),intent(in) :: coord
      integer,dimension(nnuc),intent(in)            :: iznuc 
      integer,intent(in)                            :: nfitt,fitcase 
      double precision,intent(in)                   :: v(nfitt)

      integer, dimension(3), intent(in) :: nr
      double precision, dimension(3), intent(in) :: ri, dr

      double precision :: rho,x,y,z
      integer :: i,j,k,rc
      double precision,allocatable,dimension(:) :: hgf

      allocate(hgf(nfitt),stat=rc)
      if (rc.ne.0) call errmsg('error allocating memeory')

      open(18,file=filename,status='unknown')

      write(18,fmt='(A4)') 'text'
      write(18,fmt='(A4)') 'text'
      write(18,fmt='(I5,3F12.6)')  nnuc, ri(1),  ri(2), ri(3)
      write(18,fmt='(I5,3F12.6)') nr(1), dr(1),   0.0 , 0.0 
      write(18,fmt='(I5,3F12.6)') nr(2),   0.0,  dr(2), 0.0 
      write(18,fmt='(I5,3F12.6)') nr(3),   0.0,    0.0, dr(3)
      do i=1,nnuc
        write(18,fmt='(I5,F12.6,3F12.6)')iznuc(i),0.0,coord(1:3,i)
      enddo

      do k=1, nr(1)
        x = ri(1) + (k-1)*dr(1)
        do j=1, nr(2)
          y = ri(2) + (j-1)*dr(2)
          do i=1, nr(3)
            z = ri(3) + (i-1)*dr(3)
            call fit_dens(rho,hgf,x,y,z,v,nfitt,fitcase)
            ! srampino: check this
            if (rho > 100.0d0) rho=100.0d0
            write(18,fmt='(ES20.12E3)') rho
          enddo
        enddo
      enddo

      deallocate(hgf)
      close(18)

      end subroutine cubefile 


      subroutine cubefile_nofit (nnuc,coord,iznuc,xmin,ymin,zmin,
     +                           nstep,step,filename)

      implicit none

      character(len=20),intent(in)                  :: filename
      integer,intent(in)                            :: nnuc,nstep  
      double precision,dimension(3,nnuc),intent(in) :: coord
      integer,dimension(nnuc),intent(in)            :: iznuc 
      double precision,intent(in)                   :: step,
     +                                                 xmin,ymin,zmin

      double precision :: grho(3),rho,x,y,z
      double precision,allocatable,dimension(:,:) :: temprho 
      double precision,dimension(3) :: mincube
      integer :: i,j,k


      mincube(1) = xmin
      mincube(2) = ymin
      mincube(3) = zmin

      open(18,file=filename,status='unknown')

      write(18,fmt='(A4)') 'text'
      write(18,fmt='(A4)') 'text'
      write(18,fmt='(I5,3F12.6)') nnuc, mincube 
      write(18,fmt='(I5,3F12.6)') nstep, step ,0 ,0 
      write(18,fmt='(I5,3F12.6)') nstep, 0,  step, 0 
      write(18,fmt='(I5,3F12.6)') nstep, 0,     0, step 
      do i=1,nnuc
        write(18,fmt='(I5,F12.6,3F12.6)')iznuc(i),iznuc(i),coord(1:3,i)
      enddo

      allocate(temprho(nstep,nstep))

      x = mincube(1) 
      do k=1,nstep
        y = mincube(2) 
        do j=1,nstep
          z = mincube(3) 
          do i=1,nstep
            call density(grho,rho,x,y,z)
            z = z + step
            if (rho > 100.0d0) rho=100.0d0;
            temprho(i,j)=rho
          enddo
          y = y + step
        enddo
        write(18,fmt='(6E13.5)') temprho(:,:)
        x = x + step
      enddo

      deallocate(temprho)
      close(18)

      end subroutine cubefile_nofit
