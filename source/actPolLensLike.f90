! CMB Lensing Likelihood Module for ACTPol
! ========================================
! actplens = ACTPolLensLike
!
! Send correspondence to: mathewm@astro.princeton.edu


! All Cls, whether for kappa or CMB, are assumed to start at l=0
! and step by 1.
!
!
!


module actPolLensLike

  ! CosmoMC modules
  use MatrixUtils
  use settings
  use CosmologyTypes
  use CosmoTheory
  use Calculator_Cosmology
  use Likelihood_Cosmology


  implicit none

  integer, parameter :: dl = KIND(1.d0)


  ! This structure will hold the theory spectra that we get from CosmoMC/CAMB
  type theoryCls

     ! convergence spectra, will need to convert from phi_phi
     real(8), dimension(:), allocatable :: Clkk
     ! lensed CMB spectra
     real(8), dimension(:), allocatable :: lClTT
     real(8), dimension(:), allocatable :: lClTE
     real(8), dimension(:), allocatable :: lClEE

  end type theoryCls

  ! A class that derives from a Cosmomc class that encapsulates all the likelihood data and methods
  type, extends(TCosmoCalcLikelihood) :: actplensLikelihood
     ! Estimated minimum variance Clkk in each bin from data
     real(8), dimension(:), allocatable :: lmins,lmaxes,ells, estClkk
     ! Fiducial Cls, unbinned 
     type(theoryCls) :: fCl


     ! Inverted bin-bin covariance matrix
     real(8), dimension(:,:), allocatable :: sigInv
     integer nBins ! number of Clkk_estimate bins
     integer ellMax ! this is the maximum cmb ell to use
     integer bintype ! type of binning, see ini file


     real(8), dimension(:,:), allocatable :: binMatrix  ! N0 correction matrix
     real(8), dimension(:,:), allocatable :: dlogRxRydC  ! N0 correction matrix
     real(8), dimension(:), allocatable :: N1Template ! N1 correction fiducial
     real(8), dimension(:), allocatable :: fixedn0 ! if fixing N0 correction

     real(8) :: normScale ! multiply estimated bandpowers and sqrt(covariance) by a number for debugging

     logical do_likemod ! master switch for corrections
     logical do_n0 ! n0 correction
     logical fix_n0 ! fixed n0 correction
     logical do_n1 ! n1 correction

   contains
     ! These override the main public functions for the likelihood
     procedure :: LogLike => actplens_getLnLike
     procedure :: actplens_init
     ! Some internal functions
     procedure, private :: actplens_calcLike
  end type actplensLikelihood

  
  contains

    ! ===== This subroutine interfaces with CosmoMC, called from DataLikelihoods.f90 =====
    subroutine ACTPolLensingLikelihood_Add(LikeList, Ini)
      class(TLikelihoodList) :: LikeList
      class(TSettingIni) :: ini

      Type(actplensLikelihood), pointer :: this
      integer ix

      ! Read from ini 

      write(*,*) 'In ACTPolLensingLikelihood_Add...'
      if (Ini%Read_Logical('use_actplens',.false.)) then
         write(*,*) 'Loading ACTPol Lensing data set...'
         allocate(this)

         this%name = "actpol_lensing_likelihood"
         this%do_likemod = Ini%Read_Logical('actplens_do_like_corrections',.true.)
         this%do_n0 = Ini%Read_Logical('actplens_do_n0_correction',.true.)
         this%do_n1 = Ini%Read_Logical('actplens_do_n1_kappa_correction',.true.)
         this%ellmax = Ini%Read_Int('lmin_store_all_cmb',3000)

         this%fix_n0 = Ini%Read_Logical('fix_n0_fiducial',.false.)
         this%normScale = Ini%Read_Double('normScale',1.d0)
     

         call this%actplens_init(Ini)
         call LikeList%Add(this)
      end if
      

      ! Add to list of likelihoods

    end subroutine ACTPolLensingLikelihood_Add


    subroutine actplens_init(this, Ini)
      ! This function will initialize the module by loading Cls and Cinv
      ! from data files

      class(actplensLikelihood) this
      class(TSettingIni) :: Ini
      integer nout, ix


      character(LEN=:), allocatable :: dataRoot




      dataRoot = trim(Ini%Read_String('actplens_data_root',NotFoundFail=.true.))//trim(Ini%Read_String('actplens_dataset_name',NotFoundFail=.true.))


      !!! Read inverse covmat
      write(*,*) 'Loading invcovmat...',dataRoot//"_siginv.txt"
         
      call readMatrix(dataRoot//"_siginv.txt",this%sigInv)
      this%nBins = size(this%sigInv,1)


      !!! Read estimated Cls

      write(*,*) 'Reading Cls'


      call readCls(dataRoot//"_cls.txt",this%lmins,this%lmaxes,this%ells,this%estClkk,nout)
      

      write (*,*) 'Reading bin matrix'
      call readRectMatrix(dataRoot//"_binmatrix.txt",this%nBins,this%binMatrix)
      
      !!! Likelihood corrections

      if (this%do_likemod) then
         write (*,*) 'Reading fiducial cls'
         call readCol(dataRoot//"_fidkk.txt",this%fCl%Clkk,nout)

         ! Load derivatives of Cls around fiducial
         if (this%do_n0) then


            if (this%fix_n0) then
               write (*,*) 'Reading bin matrix'
      
               call readCol(dataRoot//"_fidn0corr_"//Ini%Read_String('fix_n0_fiducial_correction_name',NotFoundFail=.true.)//".txt",this%fixedn0,nout)

            else
      
               call readCol(dataRoot//"_fidtt.txt",this%fCl%lClTT,nout)
               call readCol(dataRoot//"_fidte.txt",this%fCl%lClTE,nout)
               call readCol(dataRoot//"_fidee.txt",this%fCl%lClEE,nout)
               
            
               !dlogRRCDir = Ini%Read_String('actplens_n0_correction_matrices_path',NotFoundFail=.true.)
               write(*,*) 'Reading N0 derivatives ' !// dlogRRCDir
               call readRectMatrix(dataRoot//"_dn0.txt",this%nBins,this%dlogRxRydC)
            end if 

         end if

         ! N1 correction

         if (this%do_n1) then



            write(*,*) 'Reading N1 template ' !// n1Path
            call readCol(dataRoot//"_n1coadd.txt",this%N1Template,nout)
         

         end if

      end if



  

    end subroutine actplens_init


    function actplens_getLnLike(this,CMB, Theory, DataParams)
      use MatrixUtils
      Class(CMBParams) CMB
      Class(TCosmoTheoryPredictions), target :: Theory
      class(actplensLikelihood) this
      real(mcp) :: DataParams(:)
      real(mcp) actplens_getLnLike

      real(8), dimension(0:this%nBins-1) :: binnedFidClkk
      real(8), dimension(0:this%nBins-1) :: binnedTheoryClkk
      real(8), dimension(0:this%nBins-1) :: binnedMVTheoryClkk
      real(8), dimension(0:this%nBins-1) :: corrN0
      real(8), dimension(0:this%nBins-1) :: corrN1
      real(8) :: predperN0
      real(8) :: perN0
      real(8) :: perN1
      real(8) :: lfact
      real(8) :: b
      
      real(8) :: omegam


      type(theoryCls) :: nCl


      integer lmax, l, ix, lfrom, lto,lfloor, n1lmin,n1lmax
      
      real(8), dimension(:), allocatable :: N0Corr
      real, dimension(:), allocatable :: tElls
      real, dimension(:), allocatable :: Clvector
      real, dimension(:), allocatable :: fClvector
      real, dimension(:), allocatable :: diffCls
      real(8), parameter :: TCMB = 2.7255e6




      lmax = int(maxval(this%lmaxes))
      allocate(nCl%clkk(0:lmax))
      allocate(tElls(0:lmax))

      n1lmin = 0
      n1lmax = lmax
      
      do l=0,lmax
         if (l<2) then
            b=0.
         else
            b=1.
         end if
         tElls(l) = l
         !nCl%clkk(l) = Theory%Cls(4,4)%Cl(l)*l*2.*pi/(l+1.)/4.
         nCl%clkk(l) = b*Theory%Cls(4,4)%Cl(l)*2.*pi/4.  
      end do

      
      allocate(nCl%lcltt(0:this%ellmax))
      allocate(nCl%lclte(0:this%ellmax))
      allocate(nCl%lclee(0:this%ellmax))
         





      do l=0,this%ellMax
         if (l<2) then
            b=0.
         else
            b=1.
         end if
         lfact = b*(2.*pi)/(l*(l+1.0)) / TCMB**2
         nCl%lclTT(l) = Theory%Cls(1,1)%Cl(l)*lfact
         nCl%lclEE(l) = Theory%Cls(2,2)%Cl(l)*lfact
         nCl%lclTE(l) = Theory%Cls(2,1)%Cl(l)*lfact

      end do




      binnedTheoryClkk = matmul(transpose(this%binMatrix),nCl%clkk)
      

      !! DEBUG BINNING
      !open(unit = 1 , file = "/astro/u/msyriac/repos/DerivGen/data/bfkk.dat")
      !write (1,*) binnedFidClkk
      !close(unit=1)



      !write (*,*) "Saved file"
      !stop

      !write (*,*) "Doing corrections"

      corrN0 = 0.

      if (this%do_n0 .and. this%do_likemod) then 
         if (this%fix_n0) then
            corrN0 = this%fixedn0
         else
            binnedFidClkk = matmul(transpose(this%binMatrix),this%fCl%clkk)

            allocate(N0Corr(0:this%nBins-1))

         
            allocate(Clvector(3*(this%ellMax+1)))
            allocate(fClvector(3*(this%ellMax+1)))
            allocate(diffCls(3*(this%ellMax+1)))


            Clvector = [nCl%lClTT,nCl%lClTE,nCl%lClEE]
            fClvector = [this%fCl%lClTT(0:this%ellmax),this%fCl%lClTE(0:this%ellmax),this%fCl%lClEE(0:this%ellmax)]
            diffCls = Clvector - fClvector
            N0corr = matmul(transpose(this%dlogRxRydC) , diffCls)
            !! DEBUG BINNING
            ! open(unit = 1 , file = "/astro/u/msyriac/repos/DerivGen/data/diffcls.dat")
            ! write (1,*) diffcls
            ! close(unit=1)
            ! write (*,*) "Saved file"
            ! open(unit = 1 , file = "/astro/u/msyriac/repos/DerivGen/data/n0save.dat")
            ! write (1,*) N0corr
            ! close(unit=1)
            ! write (*,*) "Saved file"            
            ! stop
            
            corrN0 = N0corr*binnedFidClkk
            !perN0 = maxval(abs((100.*((corrN0))/(binnedTheoryClkk))))
            !write (*,*) "N0 % correction ", perN0
         end if
      end if
      
      


      corrN1 = 0.
      if (this%do_n1  .and. this%do_likemod) then

         corrN1 = this%N1Template * (((sum(tElls(n1lmin:n1lmax)*nCl%clkk(n1lmin:n1lmax)))/(sum(tElls(n1lmin:n1lmax)*this%fCl%clkk(n1lmin:n1lmax))))-1.)
         !perN1 = maxval(abs(100.*corrN1/binnedTheoryClkk))
         !write (*,*) "N1 % correction ", perN1
      end if
      
      binnedMVTheoryClkk = (binnedTheoryClkk + corrN0 + corrN1)
      


      
      
      
      ! ! return to CosmoMC the calculated -2LnLike
      actplens_getLnLike = this%actplens_calcLike(binnedMVTheoryClkk)

    end function actplens_getLnLike



    real(mcp) function actplens_calcLike(this, binnedMVTheoryClkk)
      class(actplensLikelihood) this


      real(8), dimension(0:this%nBins-1) :: binnedMVTheoryClkk
      real(8), dimension(0:this%nBins-1) :: diff
      real(8), dimension(0:this%nBins-1) :: interm


      ! Calculate and return Like

      ! This just implements 
      ! -2LnLike = (data - model) * Cinv * (data - model)^T
      diff = binnedMVTheoryClkk-(this%estClkk*this%normScale)
      actplens_calcLike = Matrix_QuadForm(this%sigInv/this%normScale/this%normScale,diff) / 2.



    end function actplens_calcLike


    !!! GENERIC UTILITIES FOR READING DATA FROM FILES


    subroutine readCls(fname,lmin,lmax,ells,Cls,nout,nmax)
      implicit none

      character (len=*) :: fname
      integer, intent(out) :: nout
      integer, optional, intent(in) :: nmax
      integer n, ix, io,nm
      integer, parameter :: funit = 100
      logical fexists
      real(8), dimension(:), allocatable :: lmin,lmax,ells, Cls
      real(8) dummy

      


      inquire(FILE=fname,EXIST=fexists)
            
      if (fexists) then

         open(unit = funit, file = fname, status = 'old', action = 'read')
         n = 0
         do
            read(funit,*,iostat=io) dummy
            if (io/=0) exit
            n = n + 1
         end do


         ! Check that user has provided nmax and if yes, that it is less than n
         nm = n
         if (present(nmax)) then
            if (nmax<n) then 
               nm = nmax
            end if
         end if
         
         nout = nm
         
         rewind(funit)

         allocate(Cls(0:nm-1))
         allocate(lmin(0:nm-1))
         allocate(lmax(0:nm-1))
         allocate(ells(0:nm-1))

         do ix = 0,nm-1
            read(funit,*) lmin(ix),lmax(ix),ells(ix),Cls(ix)
         end do


         close(funit)

      else
         write (*,*) "ERROR: Could not find file ", fname
         write (*,*) "Exiting."
         stop
         
      end if

    end subroutine readCls


    subroutine readClsSingle(fname,ells,Cls,nout,nmax)
      implicit none

      character (len=*) :: fname
      integer, intent(out) :: nout
      integer, optional, intent(in) :: nmax
      integer n, ix, io,nm
      integer, parameter :: funit = 100
      logical fexists
      real(8), dimension(:), allocatable :: ells, Cls
      real(8) dummy

      


      inquire(FILE=fname,EXIST=fexists)
            
      if (fexists) then

         open(unit = funit, file = fname, status = 'old', action = 'read')
         n = 0
         do
            read(funit,*,iostat=io) dummy
            if (io/=0) exit
            n = n + 1
         end do


         ! Check that user has provided nmax and if yes, that it is less than n
         nm = n
         if (present(nmax)) then
            if (nmax<n) then 
               nm = nmax
            end if
         end if
         
         nout = nm
         
         rewind(funit)

         allocate(Cls(0:nm-1))
         allocate(ells(0:nm-1))

         do ix = 0,nm-1
            read(funit,*) ells(ix),Cls(ix)
         end do


         close(funit)

      else
         write (*,*) "ERROR: Could not find file ", fname
         write (*,*) "Exiting."
         stop
         
      end if

    end subroutine readClsSingle


    subroutine readCol(fname,Cls,nout,nmax)
      implicit none

      character (len=*) :: fname
      integer, intent(out) :: nout
      integer, optional, intent(in) :: nmax
      integer n, ix, io,nm
      integer, parameter :: funit = 100
      logical fexists
      real(8), dimension(:), allocatable :: Cls
      real(8) dummy

      


      inquire(FILE=fname,EXIST=fexists)
            
      if (fexists) then

         open(unit = funit, file = fname, status = 'old', action = 'read')
         n = 0
         do
            read(funit,*,iostat=io) dummy
            if (io/=0) exit
            n = n + 1
         end do


         ! Check that user has provided nmax and if yes, that it is less than n
         nm = n
         if (present(nmax)) then
            if (nmax<n) then 
               nm = nmax
            end if
         end if
         
         nout = nm
         
         rewind(funit)

         allocate(Cls(0:nm-1))

         do ix = 0,nm-1
            read(funit,*) Cls(ix)
         end do


         close(funit)

      else
         write (*,*) "ERROR: Could not find file ", fname
         write (*,*) "Exiting."
         stop
         
      end if

    end subroutine readCol


    subroutine readMatrix(fname,mat)

      !
      ! Warning: This function will assume the number of columns
      ! is equal to the number of rows in the file.
      !

      implicit none

      character (len=*) :: fname
      integer n, ix, io
      integer, parameter :: funit = 100
      logical fexists
      real(8), dimension(:,:), allocatable :: mat
      real dummy


      inquire(FILE=fname,EXIST=fexists)
            
      if (fexists) then

         open(unit = funit, file = fname, status = 'old', action = 'read')
         n = 0
         do
            read(funit,*,iostat=io) dummy
            if (io/=0) exit
            n = n + 1
         end do

         rewind(funit)

         allocate(mat(0:n-1,0:n-1))


         do ix = 0,n-1
            read(funit,*) mat(ix,:)
         end do


         close(funit)

      else
         write (*) "ERROR: Could not find file", fname
         write (*) "Exiting."
         stop
         
      end if

    end subroutine readMatrix

    subroutine readRectMatrix(fname,nbins,mat)

      !
      ! Warning: This function takes as input number of columns, and determines
      ! number of rows
      !

      implicit none

      integer nbins
      character (len=*) :: fname
      integer n, ix, io
      integer, parameter :: funit = 100
      logical fexists
      real(8), dimension(:,:), allocatable :: mat
      real dummy


      inquire(FILE=fname,EXIST=fexists)
            
      if (fexists) then

         open(unit = funit, file = fname, status = 'old', action = 'read')
         n = 0
         do
            read(funit,*,iostat=io) dummy
            if (io/=0) exit
            n = n + 1
         end do

         rewind(funit)

         allocate(mat(0:n-1,0:nbins-1))


         do ix = 0,n-1
            read(funit,*) mat(ix,:)
            !write(*,*) mat(ix,:)
         end do


         close(funit)

      else
         write (*) "ERROR: Could not find file", fname
         write (*) "Exiting."
         stop
         
      end if

    end subroutine readRectMatrix


end module actPolLensLike



    !        MMM.           .MMM
    !        MMMMMMMMMMMMMMMMMMM
    !        MMMMMMMMMMMMMMMMMMM      ___________________________________
    !       MMMMMMMMMMMMMMMMMMMMM    |                                   |
    !      MMMMMMMMMMMMMMMMMMMMMMM   | Octokitty is judging your code.   |
    !     MMMMMMMMMMMMMMMMMMMMMMMM   |_   _______________________________|
    !     MMMM::- -:::::::- -::MMMM    |/
    !      MM~:~   ~:::::~   ~:~MM
    ! .. MMMMM::. .:::+:::. .::MMMMM ..
    !       .MM::::: ._. :::::MM.
    !          MMMM;:::::;MMMM
    !   -MM        MMMMMMM
    !   ^  M+     MMMMMMMMM
    !       MMMMMMM MM MM MM
    !            MM MM MM MM
    !            MM MM MM MM
    !         .~~MM~MM~MM~MM~~.
    !      ~~~~MM:~MM~~~MM~:MM~~~~
    !     ~~~~~~==~==~~~==~==~~~~~~
    !      ~~~~~~==~==~==~==~~~~~~
    !          :~==~==~==~==~~
