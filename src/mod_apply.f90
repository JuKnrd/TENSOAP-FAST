module apply

 real*8, allocatable :: xyz(:,:,:),cell(:,:,:)
 character(len=4), allocatable :: atname(:,:)
 character(len=1000), allocatable :: comment(:)
 integer, allocatable :: natoms(:)
 integer nframes,nlines,natmax

 contains

!****************************************************************************************************************
subroutine get_xyz(fname,periodic)
 implicit none

  character(len=100) fname
  character(len=1000) c1
  integer ios,i,j,frnum,k,ii
  logical periodic

    ! Read in XYZ file
    open(unit=31,file=fname,status='old')
    nlines = 0
    do
     read(31,*,iostat=ios)
     if (ios.ne.0) exit
     nlines = nlines + 1
    enddo
    close(31)

    open(unit=31,file=fname,status='old')
    nframes = 0
    natmax = 0
    do i=1,nlines
     read(31,*,iostat=ios) frnum
     if (ios.eq.0) then
      nframes = nframes + 1
      natmax = max(natmax,frnum)
     endif
    enddo
    close(31)

    allocate(xyz(nframes,natmax,3),atname(nframes,natmax),natoms(nframes),comment(nframes))
    xyz(:,:,:) = 0.d0
    atname(:,:) = ''
    open(unit=31,file=fname,status='old')
    k = 0
    do i=1,nframes
     read(31,*,iostat=ios) frnum
     if (ios.eq.0) then
      k = k + 1
      natoms(k) = frnum
      read(31,'(A)') comment(k)
      do j=1,natoms(k)
       read(31,*) atname(k,j),(xyz(k,j,ii),ii=1,3)
      enddo
     else
      write(*,*) 'At frame ',k,':'
      stop 'ERROR: there should be an atom number here!'
     endif
    enddo

    ! Get cell data
    allocate(cell(nframes,3,3))
    if (.not.periodic) then
     cell(:,:,:) = 0.d0
    else
     do i=1,nframes
      ios = index(comment(i),'Lattice')
      if (ios.eq.0) stop 'ERROR: input file is not periodic!'
      c1 = comment(i)
      c1 = c1(ios:len(c1))
      ios = index(c1,'"')
      c1 = c1(ios+1:len(c1))
      ios = index(c1,'"')
      c1 = c1(1:ios-1)
      read(c1,*) cell(i,1,1),cell(i,1,2),cell(i,1,3),cell(i,2,1),cell(i,2,2),cell(i,2,3),cell(i,3,1),cell(i,3,2),cell(i,3,3)
     enddo
    endif

    return

end subroutine
!****************************************************************************************************************

end module
