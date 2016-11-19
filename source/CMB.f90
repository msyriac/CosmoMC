
    module CMBLikelihoods
    use Likelihood_Cosmology
    use CosmologyTypes
    use CosmoTheory
    use CMBLikes
    implicit none
    private

    type, extends(TCMBLikelihood) :: TCMBSZLikelihood
     ! (inherited) tag set from "cmb_dataset[tag] =" in input file
        real(mcp), pointer, dimension(:) :: sz_template
    contains
    procedure :: ReadSZTemplate
    procedure :: ReadParams => TCMBSZLikelihood_ReadParams
    end type TCMBSZLikelihood

#ifdef WMAP
    type, extends(TCMBSZLikelihood) :: TWMAPLikelihood
    contains
    procedure :: ReadParams => TWMAPLikelihood_ReadParams
    procedure :: LogLike => TWMAPLikelihood_LogLike
    end type TWMAPLikelihood
#endif

!Erminia:add ACTPol

 type, extends(TCMBLikelihood) :: CMBAPLikelihood
    contains
    procedure :: ReadParams => CMBAPLikelihood_ReadParams
    end type CMBAPLikelihood

#ifdef ACTPol
    type, extends(CMBAPLikelihood) :: ACTPolLikelihood
    contains
    procedure :: ReadParams => ACTPolLikelihood_ReadParams
    procedure :: LogLike => ACTPolLikelihood_LogLike
    end type ACTPolLikelihood
#endif
!!

!Erminia:add ACT/SPT

 type, extends(TCMBLikelihood) :: CMBACTSPTLikelihood
    contains
    procedure :: ReadParams => CMBACTSPTLikelihood_ReadParams
    end type CMBACTSPTLikelihood

#ifdef ACTSPT
    type, extends(CMBACTSPTLikelihood) :: ACTSPTLikelihood
    contains
    procedure :: ReadParams => ACTSPTLikelihood_ReadParams
    procedure :: LogLike => ACTSPTLikelihood_LogLike
    end type ACTSPTLikelihood
#endif
!!

    public CMBLikelihood_Add
    contains


    subroutine CMBLikelihood_Add(LikeList, Ini)
#ifdef CLIK
    use cliklike
#endif
#ifdef NONCLIK
    use noncliklike
#endif
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    class(TCMBLikelihood), pointer  :: like
    integer  i
    Type(TSettingIni) :: DataSets

    call Ini%TagValuesForName('cmb_dataset', DataSets, filename=.true.)

    do i= 1, DataSets%Count
        if (DataSets%Name(i) == 'WMAP') then
#ifdef WMAP
            allocate(TWMAPLikelihood::like)
            like%name = Datasets%Value(i)
            select type(like)
            class is (TWMAPLikelihood)
            end select
#else
            call MpiStop('Set WMAP directory in Makefile to compile with WMAP')
#endif
!! Erminia
!!
        else if (DataSets%Name(i) == 'ACTPol') then
#ifdef ACTPol 
            allocate(ACTPolLikelihood::like)
            like%name = Datasets%Value(i)
            select type(like)
            class is (ACTPolLikelihood)
            like%Tag = DataSets%Name(i)
            end select
#else
            call MpiStop('Set ACTPol directory in Makefile to compile with ACTPol')
#endif
!!

!! Erminia
!!
        else if (DataSets%Name(i) == 'ACTSPT') then
#ifdef ACTSPT
            allocate(ACTSPTLikelihood::like)
            like%name = Datasets%Value(i)
            select type(like)
            class is (ACTSPTLikelihood)
            like%Tag = DataSets%Name(i)
            end select
#else
            call MpiStop('Set ACTSPT directory in Makefile to compile with ACT/SPT')
#endif
!!

        else
            allocate(TCMBLikes::like)
            call like%ReadDatasetFile(Datasets%Value(i))
        end if
        like%Tag = DataSets%Name(i)
        call like%ReadParams(Ini)
        call Ini%Read(Ini%NamedKey('cmb_dataset_speed',DataSets%Name(i)),like%speed)

        call LikeList%Add(like)
    end do
    if (Feedback > 1 .and. DataSets%Count>0 ) write (*,*) 'read CMB data sets'

#ifdef CLIK
    Use_clik = Ini%Read_Logical('use_clik',.false.)
    if (Use_clik) then
        call clik_readParams(LikeList, Ini)
    end if
#else
    if (Ini%Read_Logical('use_clik',.false.)) call MpiStop('compile with CLIK to use clik - see Makefile')
#endif
#ifdef NONCLIK
    call nonclik_readParams(LikeList, Ini)
#endif

    end subroutine CMBLikelihood_Add


    subroutine TCMBSZLikelihood_ReadParams(this, Ini)
    class(TCMBSZLikelihood) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: SZTemplate
    real(mcp) :: SZscale = 1

    SZTemplate = Ini%Read_String(Ini%NamedKey('cmb_dataset_SZ',this%tag))
    if (SZTemplate/='') then
        call Ini%Read(Ini%NamedKey('cmb_dataset_SZ_scale',this%tag),SZScale)
        call this%ReadSZTemplate(SZTemplate,SZScale)
        call this%loadParamNames(trim(DataDir)//'WMAP.paramnames')
    end if

    call this%TCMBLikelihood%ReadParams(Ini)

    end subroutine TCMBSZLikelihood_ReadParams

    subroutine ReadSZTemplate(this, aname, ascale)
    class(TCMBSZLikelihood) :: this
    character(LEN=*), intent(IN) :: aname
    real(mcp), intent(in) :: ascale
    integer l, status
    real(mcp) sz
    Type(TTextFile) :: F

    allocate(this%sz_template(2:this%cl_lmax(CL_T,CL_T)))
    this%sz_template = 0
    call F%Open(aname)
    do
        read(F%unit,*,iostat=status) l, sz
        if (status/=0) exit
        if (l>=2 .and. l<=this%cl_lmax(CL_T,CL_T)) this%sz_template(l) = ascale * sz
    end do
    call F%Close()

    end subroutine ReadSZTemplate

#ifdef WMAP
    subroutine TWMAPLikelihood_ReadParams(this, Ini)
    use WMAP_OPTIONS
    class(TWMAPLikelihood) :: this
    class(TSettingIni) :: ini

    use_TT_beam_ptsrc = Ini%read_Logical('use_WMAP_TT_beam_ptsrc', .true.)
    use_TE = Ini%read_Logical('use_WMAP_TE',.true.)
    use_TT = Ini%read_Logical('use_WMAP_TT',.true.)
!Erminia
    use_lowl_TT = Ini%read_Logical('use_WMAP_lowl_TT',.true.)
    use_lowl_pol = Ini%read_Logical('use_WMAP_lowl_pol',.true.)
!
    if (MPIRank==0) print *, 'WMAP options (beam TE TT)', use_TT_beam_ptsrc, use_TE, use_TT
    allocate(this%cl_lmax(CL_B,CL_B), source=0)
    this%cl_lmax(CL_T,CL_T) = ttmax
    this%cl_lmax(CL_E,CL_T) = temax
    this%cl_lmax(CL_E,CL_E) = max(gibbs_ell_max,lowl_max)
    this%cl_lmax(CL_B,CL_B) = max(gibbs_ell_max,lowl_max)

    call this%TCMBSZLikelihood%ReadParams(Ini)

    end subroutine TWMAPLikelihood_ReadParams

    function TWMAPLikelihood_LogLike(this, CMB, Theory, DataParams) result(logLike)
    use wmap_likelihood_9yr
    use WMAP_OPTIONS
    use WMAP_UTIL
    Class(TWMAPLikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) logLike
    real(8) :: likes(num_WMAP),like_tot
    real(mcp) CLTT(2:this%cl_lmax(1,1))

    CLTT = Theory%Cls(1,1)%CL(2:this%cl_lmax(1,1)) + DataParams(1)*this%sz_template
    likes=0
    call wmap_likelihood_compute(CLTT,Theory%Cls(2,1)%CL(2:),Theory%Cls(2,2)%CL(2:),Theory%Cls(3,3)%CL(2:),likes)
    !call wmap_likelihood_error_report

    if (wmap_likelihood_ok) then
        LogLike = sum(likes)
    else
        LogLike = LogZero
    endif
    end function TWMAPLikelihood_LogLike
#endif

!! Erminia 

    subroutine CMBAPLikelihood_ReadParams(this, Ini)
    class(CMBAPLikelihood) :: this
    class(TSettingIni) :: Ini
    call this%loadParamNames(trim(DataDir)//'ACTPol.paramnames')
    call this%TCMBLikelihood%ReadParams(Ini)
    end subroutine CMBAPLikelihood_ReadParams

#ifdef ACTPol
    subroutine ACTPolLikelihood_ReadParams(this, Ini)
    use ACTPol_CMBonly
    class(ACTPolLikelihood) :: this
    class(TSettingIni) :: ini
    allocate(this%cl_lmax(CL_B,CL_B), source=0)
    this%cl_lmax(CL_T,CL_T) = tt_lmax
    this%cl_lmax(CL_E,CL_T) = tt_lmax
    this%cl_lmax(CL_E,CL_E) = tt_lmax
    call this%CMBAPLikelihood%ReadParams(Ini)
    end subroutine ACTPolLikelihood_ReadParams

  function ACTPolLikelihood_LogLike(this, CMB, Theory, DataParams) result(logLike)
    use ACTPol_CMBonly
    Class(ACTPolLikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) logLike
    real(mcp) like
    logical :: init_actpol = .true.
    real(mcp), dimension(1) :: fgp

    if (init_actpol) then
      call actpol_like_init
      if (Feedback>0) write(*,*) 'reading ACTPol data'
      init_actpol = .false.
    end if

    fgp(:)=DataParams(:)

    like=0.d0

    call actpol_calc_like(like,Theory%Cls(1,1)%CL(2:),Theory%Cls(2,1)%CL(2:),Theory%Cls(2,2)%CL(2:),fgp(1))

    LogLike = like
    write(*,*) 'ACTPolLnLike=', LogLike
   end function ACTPolLikelihood_LogLike
#endif

!! Erminia 

    subroutine CMBACTSPTLikelihood_ReadParams(this, Ini)
    class(CMBACTSPTLikelihood) :: this
    class(TSettingIni) :: Ini
    call this%loadParamNames(trim(DataDir)//'ACTSPT.paramnames')
    call this%TCMBLikelihood%ReadParams(Ini)
    end subroutine CMBACTSPTLikelihood_ReadParams

#ifdef ACTSPT
    subroutine ACTSPTLikelihood_ReadParams(this, Ini)
    use actlite_3yr_like
    class(ACTSPTLikelihood) :: this
    class(TSettingIni) :: ini
    allocate(this%cl_lmax(CL_B,CL_B), source=0)
    this%cl_lmax(CL_T,CL_T) = tt_lmax
    call this%CMBACTSPTLikelihood%ReadParams(Ini)
    end subroutine ACTSPTLikelihood_ReadParams

  function ACTSPTLikelihood_LogLike(this, CMB, Theory, DataParams) result(logLike)
    use actlite_3yr_like
    Class(ACTSPTLikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) logLike
    real(8) :: like
    logical :: init_actspt = .true.

    if (init_actspt) then
      call act_likelihood_init
      if (Feedback>0) write(*,*) 'reading ACT/SPT data'
      init_actspt = .false.
    end if

    like=0.d0

    call act_likelihood_compute(Theory%Cls(1,1)%CL(2:),like)
    LogLike = like
    write(*,*) 'ACTSPTLnLike=', LogLike
   end function ACTSPTLikelihood_LogLike
#endif

    end module CMBLikelihoods
