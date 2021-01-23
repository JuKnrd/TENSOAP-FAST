program sagpr_apply
    use apply
    use io
    USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
    implicit none

    type(SAGPR_Model) :: GPR
    character(len=100), allocatable :: arg(:),keys(:)
    integer i,j,k,l,ios,numkeys,port,socket,inet,cbuf,nargs
    integer nat
    character(len=100) ofile,model,fname,sock_arg(3)
    character(len=1024) hostname
    character(len=12) header
    character(len=2048) initbuffer
    integer, parameter :: MSGLEN=12
    integer t1,t2,cr,ts,tf
    real*8 rate,mtxbuf(9),cell_h(3,3),virial(3,3)
    real*8, allocatable :: msgbuffer(:)
    real*8, parameter :: bohrtoangstrom = 0.52917721d0
    logical readnext,use_socket,hasdata

!************************************************************************************
! GET COMMAND-LINE ARGUMENTS
!************************************************************************************

    ! Get input arguments
    nargs = iargc()
    allocate(arg(nargs+1))
    do i=1,nargs
     call getarg(i,arg(i))
    enddo
    arg(nargs+1) = 'NULL'

    ! Parse these arguments
    numkeys = 8
    allocate(keys(numkeys))
    model = ''
    ofile = 'prediction.out'
    fname = ''
    use_socket = .false.
    inet = 1
    GPR%verbose = .false.
    GPR%atomic = .false.
    keys = (/'-m  ','-o  ','-f  ','-s  ','-u  ','-v  ','-a  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-m') read(arg(i+1),'(A)') model
     if (arg(i).eq.'-o') read(arg(i+1),'(A)') ofile
     if (arg(i).eq.'-f') read(arg(i+1),'(A)') fname
     if (arg(i).eq.'-v') GPR%verbose=.true.
     if (arg(i).eq.'-a') GPR%atomic=.true.
     if (arg(i).eq.'-s') then
      use_socket = .true.
      readnext = .true.
      l = 0
      do k=i+1,nargs
       do j=1,numkeys
        readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keys(j))))
       enddo
       if (readnext) then
        l = l + 1
        read(arg(k),*) sock_arg(l)
       endif
      enddo
     endif
     if (arg(i).eq.'-u') inet = 0
    enddo
    deallocate(arg)

    ! Check for arguments that are required
    if (model.eq.'') stop 'ERROR: model file required!'
    if (fname.eq.'') stop 'ERROR: file name required!'

!************************************************************************************
! GET MODEL
!************************************************************************************

    call get_model(GPR,model)

!************************************************************************************
! READ IN DATA AND MAKE PREDICTIONS
!************************************************************************************

    ! Open xyz file
    open(unit=15,file=trim(adjustl(fname)),status="old",access="stream",form="formatted",iostat=ios)
    if (ios.ne.0) stop 'ERROR: input file does not exist!'
    if (use_socket) then
     ! Get atom names
     call read_frame(GPR%frames,15,GPR%verbose,GPR%periodic)
     ! Set up socket
     hostname = trim(adjustl(sock_arg(1)))//achar(0)
     read(sock_arg(2),*) port
     call open_socket(socket,inet,port,hostname)
    endif

    ! Initialize the system clock
    call system_clock(count_rate=cr)
    rate = real(cr)

    ios = 1
    hasdata = .false.
    ! Open output file
    if (.not. use_socket) open(unit=33,file=ofile,access='stream',form='formatted')
    if (use_socket .and. GPR%atomic) open(unit=33,file=ofile,access='stream',form='formatted')
    ! Keep going through input
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      if (.not. use_socket) then
       ! Read the next frame from xyz file
       call read_frame(GPR%frames,15,GPR%verbose,GPR%periodic)

       call system_clock(t1)

       ! Do prediction
       call predict_frame(GPR,rate)

       ! Rescale predictions about the mean
       call rescale_predictions(GPR)

       ! Print predictions
       call print_predictions(GPR,33)

       call system_clock(t2)
       if (GPR%verbose) write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
       if (GPR%verbose) write(*,*)

      else
       ! Read from the socket
       call readbuffer(socket,header,MSGLEN)
       if (trim(header)=='STATUS') then
        ! The wrapper is inquiring what we are doing
        if (hasdata) then
         call writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
        else
         call writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
        endif
       elseif (trim(header)=='INIT') then
        ! The wrapper is sending an initialization string (which we will ignore)
        call readbuffer(socket, cbuf)
        call readbuffer(socket, initbuffer, cbuf)
       elseif (trim(header)=='POSDATA') then
        ! The wrapper is sending atom positions
        if (allocated(GPR%frames%xyz)) deallocate(GPR%frames%xyz,GPR%frames%natoms,GPR%frames%cell)
        call readbuffer(socket, mtxbuf, 9)  ! Cell matrix
        cell_h = reshape(mtxbuf, (/3,3/))
        cell_h = transpose(cell_h)
        call readbuffer(socket, mtxbuf, 9)
        call readbuffer(socket, cbuf)
        nat = cbuf
        GPR%frames%natmax = nat
        GPR%frames%nframes = 1
        allocate(GPR%frames%xyz(GPR%frames%nframes,GPR%frames%natmax,3),GPR%frames%natoms(GPR%frames%nframes), &
     &     msgbuffer(3*nat),GPR%frames%cell(GPR%frames%nframes,3,3))
        GPR%frames%natoms(1) = nat
        call readbuffer(socket, msgbuffer, nat*3) 
        ! Read in atoms
        do i=1,nat
         GPR%frames%xyz(1,i,:) = msgbuffer(3*(i-1)+1:3*i) * bohrtoangstrom
        enddo
        GPR%frames%cell(1,:,:) = cell_h(:,:) * bohrtoangstrom

        ! With all of the necessary data read in, do the prediction
        call system_clock(t1)

        ! Do prediction
        call predict_frame(GPR,rate)

        ! Rescale committee predictions about their mean
        call rescale_predictions(GPR)

        call system_clock(t2)
        if (GPR%verbose) write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
        if (GPR%verbose) write(*,*)

        ! We have the data
        hasdata = .true.
        deallocate(GPR%frames%xyz,GPR%frames%natoms,msgbuffer,GPR%frames%cell)
       elseif (trim(header)=='GETFORCE') then
        ! Now we send all of the information back to the wrapper
        ! We start with a lot of zeros (there is no contribution to the energy
        ! or force)
        allocate(msgbuffer(3*nat))
        msgbuffer(:) = 0.d0
        virial(:,:) = 0.d0
        call writebuffer(socket,"FORCEREADY  ",MSGLEN)
        call writebuffer(socket,0.d0)
        call writebuffer(socket,nat)
        call writebuffer(socket,msgbuffer,3*nat)
        call writebuffer(socket,reshape(virial,(/9/)),9)
        ! Send prediction; we will send only the prediction for the entire
        ! frame, rather than for each atom (the latter will be stored in an
        ! output file)
        initbuffer = " "
        write(initbuffer,*) ((sum(GPR%prediction_lm_c(:,j,k)),j=1,GPR%degen),k=1,GPR%nw)
        cbuf = len_trim(initbuffer)
        call writebuffer(socket,cbuf)
        call writebuffer(socket,initbuffer,cbuf)
        ! If we are doing atomic predictions, also print them to a file
        if (GPR%atomic) then
         write(33,*) size(GPR%prediction_lm_c,1)
         write(33,*) '# Total',((sum(GPR%prediction_lm_c(:,j,k)),j=1,GPR%degen),k=1,GPR%nw)
         do l=1,size(GPR%prediction_lm_c,1)
          write(33,*)  GPR%atname_at(l),((GPR%prediction_lm_c(l,j,k),j=1,GPR%degen),k=1,GPR%nw)
         enddo
         flush(33)
        endif
        deallocate(msgbuffer)
        hasdata = .false.
       else
        write(*,*) 'ERROR: unexpected header ',header
        stop
       endif
      endif
     endif
    enddo
    close(33)

    ! Array deallocation
    if (allocated(GPR%frames%xyz)) deallocate(GPR%frames%xyz,GPR%frames%atname, &
     &     GPR%frames%natoms,GPR%frames%comment,GPR%sparsification, &
     &     GPR%frames%cell,GPR%PS_tr_lam,GPR%wt,arg,keys,GPR%prediction_lm,GPR%PS)
    if (allocated(GPR%natoms_tr)) deallocate(GPR%natoms_tr)
    if (allocated(GPR%PS0)) deallocate(GPR%PS0)
    if (allocated(GPR%ker)) deallocate(GPR%ker)
    if (allocated(GPR%ker_lm)) deallocate(GPR%ker_lm)
    if (allocated(GPR%sparsification0)) deallocate(GPR%sparsification0,GPR%PS_tr_0)
    if (allocated(GPR%components)) deallocate(GPR%components)
    if (allocated(GPR%w3j)) deallocate(GPR%w3j)
    if (allocated(GPR%meanval_c)) deallocate(GPR%meanval_c,GPR%wt_c,GPR%prediction_lm_c)
    if (allocated(GPR%atname_at)) deallocate(GPR%atname_at)

end program
