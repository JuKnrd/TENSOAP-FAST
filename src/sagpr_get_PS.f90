program sagpr_predict
    use sagpr
    implicit none

    integer nargs,numkeys
    character(len=100), allocatable :: arg(:),keylist(:)
    integer lm,i,j,nmol,reals,bytes,degen,nenv,mu,nu,k
    integer nmax,lmax,ncut
    real*8 rcut,sg,rs(3)
    character(len=100) weights,ofile,kerfile,cen,spec,sparse,fname
    real*8, allocatable, target :: all_wt(:),raw_ker(:)
    real*8, allocatable :: wt(:),prediction_lm(:,:),ker(:,:,:,:)
    real*8 meanval

    ! Get input arguments
    nargs = iargc()
    allocate(arg(nargs+1))
    do i=1,nargs
     call getarg(i,arg(i))
    enddo
    arg(nargs+1) = 'NULL'

    ! Parse these arguments
    lm = 0
    nmax = 8
    lmax = 6
    rcut = 4.d0
    sg = 0.3d0
    cen = ''
    spec = ''
    ncut = -1
    sparse = ''
    fname = ''
    rs = (/0.d0,0.d0,0.d0/)
    ofile = ''
    keylist = (/'-lm ','-n  ','-l  ','-rc ','-sg ','-c  ','-s  ','-nc ','-sf ','-f  ','-rs ','-o  ','NULL'/)

    stop

    lm = 0
    weights = ''
    kerfile = ''
    ofile = 'prediction.out'
    numkeys = 5
    allocate(keylist(numkeys))
    keylist = (/'-lm ','-w  ','-k  ','-o  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-lm') read(arg(i+1),*) lm
     if (arg(i).eq.'-w') read(arg(i+1),*) weights
     if (arg(i).eq.'-k') read(arg(i+1),*) kerfile
     if (arg(i).eq.'-o') read(arg(i+1),*) ofile
    enddo

    ! Check for arguments that are required
    if (weights.eq.'') stop 'ERROR: weights file required!'
    if (kerfile.eq.'') stop 'ERROR: kernel file required!'

    ! Read in weights
    degen = 2*lm + 1
    open(unit=31,file=weights,status='old',access='stream',form='unformatted')
    inquire(unit=31,size=bytes)
    reals = bytes/8
    allocate(all_wt(reals),wt(reals-1))
    read(31,pos=1) all_wt
    meanval = all_wt(1)
    do i=2,reals
     wt(i-1) = all_wt(i)
    enddo
    nenv = (reals-1) / degen
    close(31)

    ! Read in kernel
    open(unit=32,file=kerfile,status='old',access='stream',form='unformatted')
    inquire(unit=32,size=bytes)
    reals = bytes/8
    allocate(raw_ker(reals))
    read(32,pos=1) raw_ker
    close(32)
    nmol = reals / (nenv * degen * degen)
    allocate(ker(nmol,nenv,degen,degen))
    k = 0
    do i=1,nmol
     do j=1,nenv
      do nu=1,degen
       do mu=1,degen
        k = k + 1
        ker(i,j,mu,nu) = raw_ker(k)
       enddo
      enddo
     enddo
    enddo

    ! Do prediction
    allocate(prediction_lm(nmol,degen))
    prediction_lm = do_prediction(ker,wt,meanval,degen,nmol,nenv)

    ! Print out predictions
    open(unit=33,file=ofile)
    do i=1,nmol
     write(33,*) (prediction_lm(i,mu),mu=1,degen)
    enddo
    close(33)

    ! Array deallocation
    deallocate(all_wt,wt,keylist,ker,prediction_lm,raw_ker)

end program
