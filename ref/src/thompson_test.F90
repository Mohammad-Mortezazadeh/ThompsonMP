program test_thompson
   USE mt19937
   USE thompson_utils
   USE machine, only: kind_phys
   USE mp_thompson
#define MPI
#ifdef _OPENMP
   USE omp_lib
#endif
use mpi_f08

   IMPLICIT NONE

   !===============================
   ! Interface variables
   integer                    :: ncol
   integer                    :: nlev
   real(4)            :: con_g, con_rd, con_eps
   logical                    :: restart
   integer                    :: imp_physics
   integer                    :: imp_physics_thompson
   ! Hydrometeors
   logical                    :: convert_dry_rho
   real(4), allocatable     :: spechum(:,:)
   real(4), allocatable     :: qc(:,:)
   real(4), allocatable     :: qr(:,:)
   real(4), allocatable     :: qi(:,:)
   real(4), allocatable     :: qs(:,:)
   real(4), allocatable     :: qg(:,:)
   real(4), allocatable     :: ni(:,:)
   real(4), allocatable     :: nr(:,:)
   ! Aerosols
   logical                    :: is_aerosol_aware
   real(4), allocatable     :: nc(:,:)
   real(4), allocatable     :: nwfa(:,:)
   real(4), allocatable     :: nifa(:,:)
   real(4), allocatable     :: nwfa2d(:)
   real(4), allocatable     :: nifa2d(:)
   ! State variables
   real(4), allocatable     :: tgrs(:,:)
   real(4), allocatable     :: prsl(:,:)
   real(4), allocatable     :: phil(:,:)
   real(4), allocatable     :: area(:)
   ! MPI information
   TYPE(MPI_Comm) :: mpicomm
   integer                    :: mpirank
   integer                    :: mpiroot
   integer                    :: mpisize
   ! Threading/blocking information
   integer                    :: threads
   ! Extended diagnostics
   logical                    :: ext_diag
   real(4), allocatable     :: diag3d(:,:,:)
   ! CCPP error handling
   character(len=256)           :: errmsg
   integer                    :: errflg
   !===============================
   ! Aerosols
   logical                     :: reset_dBZ
   logical                     :: aero_ind_fdb
   ! State variables and timestep information
   real(4), allocatable            :: omega(:,:)
   real(4)            :: dtp
   logical                    :: first_time_step
   integer                    :: istep, nsteps
   real(4)            :: dt_inner
   ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
   real(4), allocatable            :: prcp(:)
   real(4), allocatable            :: rain(:)
   real(4), allocatable            :: graupel(:)
   real(4), allocatable            :: ice(:)
   real(4), allocatable            :: snow(:)
   real(4), allocatable            :: sr(:)
   ! Radar reflectivity
   real(4), allocatable            :: refl_10cm(:,:)
   ! State variables and timestep information
   real(4), allocatable     :: phii(:,:)
   logical                    :: do_radar_ref
   logical                       :: sedi_semi
   integer                       :: decfl
   ! MPI and block information
   integer                       :: blkno
   ! SPP
   integer                   :: spp_mp
   integer                   :: n_var_spp
   real(4),           allocatable :: spp_wts_mp(:,:)
   real(4),           allocatable :: spp_prt_list(:)
   character(len=10),          allocatable :: spp_var_list(:)
   real(4),           allocatable :: spp_stddev_cutoff(:)
   ! Extended diagnostic output
   logical                  :: reset_diag3d

   !===============================
   integer, parameter :: ext_ndiag3d = 37
   !===============================

   integer :: count_rate, count_start, count_end
   real :: elapsed

   integer :: alloc_stat
   integer :: n_omp_threads, s, e, tid
   integer :: N_GPUS, gpuid
   integer, parameter :: DTEND_DIM = 12

   integer ierror
   integer :: ncol_per_thread, ncol_rem
   character(64) :: str
   !===============================
   !ccpp new variables, MOHAMMAD
   logical :: is_initialized
   logical :: merra2_aerosol_aware, fullradar_diag
   real(4), allocatable :: aerfld(:,:,:)
   integer, allocatable :: islmsk(:)
   real(4), allocatable :: max_hail_diam_sfc(:)
   logical :: cplchm
   real(4), allocatable :: pfi_lsan(:,:)
   real(4), allocatable :: pfl_lsan(:,:)


#ifdef MPI
   mpicomm = MPI_COMM_WORLD
   mpiroot = 0
   CALL MPI_INIT(ierror)
   CALL MPI_COMM_SIZE(mpicomm, mpisize, ierror)
   CALL MPI_COMM_RANK(mpicomm, mpirank, ierror)

   gpuid = mpirank
#else
   gpuid = 0
#endif
   PRINT*, 'MPI rank', mpirank

#ifdef _OPENMP
!$omp parallel
!$omp single
   n_omp_threads = omp_get_num_threads()
!$omp end single
!$omp end parallel
#else
   n_omp_threads = 1
#endif

   N_GPUS = 0

#ifdef _OPENMP
   WRITE(6,'(" Using ",i3," threads")') n_omp_threads
#endif

   !===============================
   if (COMMAND_ARGUMENT_COUNT().GE.1) THEN
      CALL GET_COMMAND_ARGUMENT(1, str)
      READ(str,*) ncol
   else
      ncol = 512
   endif
   if (COMMAND_ARGUMENT_COUNT().GE.2) THEN
      CALL GET_COMMAND_ARGUMENT(2, str)
      READ(str,*) nlev
   else
      nlev = 64
   endif
   print*, ncol, nlev
   !===============================
   !===============================
   threads = n_omp_threads
   ext_diag = .TRUE.
   is_aerosol_aware = .FALSE.
   merra2_aerosol_aware = .FALSE.!MOHAMMAD
   is_initialized = .False.!MOHAMMAD
   fullradar_diag = .False.!MOHAMMAD
   cplchm = .False.!MOHAMMAD
   convert_dry_rho = .FALSE.
   con_g = 9.80665
   con_rd = 287.05
   con_eps = 0.6219993
   restart = .FALSE.
   imp_physics = 1
   imp_physics_thompson = 1
   !===============================
   reset_dBZ = .FALSE.
   aero_ind_fdb = .FALSE.
   first_time_step = .TRUE.
   dtp = 0.01
   istep = 1
   nsteps = 20
   dt_inner = 0.02
   do_radar_ref = .TRUE.
   sedi_semi = .FALSE.
   decfl = 1
   blkno = 1
   spp_mp = 1
   n_var_spp = 2
   reset_diag3d = .TRUE.
   !===============================

#ifdef MPI
   CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

   PRINT*, "Allocating arrays"
   ALLOCATE(                        &
       spechum(ncol,nlev),          &
       qc(ncol,nlev),               &
       qr(ncol,nlev),               &
       qi(ncol,nlev),               &
       qs(ncol,nlev),               &
       qg(ncol,nlev),               &
       ni(ncol,nlev),               &
       nr(ncol,nlev),               &
       nc(ncol,nlev),               &
       nwfa(ncol,nlev),             &
       nifa(ncol,nlev),             &
       nwfa2d(ncol),                &
       nifa2d(ncol),                &
       tgrs(ncol,nlev),             &
       prsl(ncol,nlev),             &
       phil(ncol,nlev),             &
       area(ncol),                  &
       diag3d(ncol,nlev,ext_ndiag3d), &
       omega(ncol,nlev),            &
       prcp(ncol),                  &
       rain(ncol),                  &
       graupel(ncol),               &
       ice(ncol),                   &
       snow(ncol),                  &
       sr(ncol),                    &
       refl_10cm(ncol,nlev),        &
       phii(ncol,nlev+1),           &
       spp_wts_mp(ncol,nlev),       &
       spp_prt_list(ncol),          &
       spp_var_list(ncol),          &
       spp_stddev_cutoff(ncol),     &
       islmsk(ncol),                &!MOHAMMAD
       aerfld(ncol,nlev,16),        &!MOHAMMAD
       max_hail_diam_sfc(ncol),     &!MOHAMMAD
       pfi_lsan(ncol, nlev),        &!MOHAMMAD
       pfl_lsan(ncol, nlev),        &!MOHAMMAD


       STAT=alloc_stat)
   IF (alloc_stat /= 0) STOP "Error allocating arrays"

   !=============================================================
   PRINT*, "Initializing arrays"
   s = 1
   e = ncol

   CALL mt19937_real2d(spechum(s:e,:))
   spechum=spechum*0.01
   CALL mt19937_real2d(qc(s:e,:))
   CALL mt19937_real2d(qr(s:e,:))
   CALL mt19937_real2d(qi(s:e,:))
   CALL mt19937_real2d(qs(s:e,:))
   CALL mt19937_real2d(qg(s:e,:))
   CALL mt19937_real2d(ni(s:e,:))
   CALL mt19937_real2d(nr(s:e,:))
   CALL mt19937_real2d(nc(s:e,:))
   CALL mt19937_real2d(nwfa(s:e,:))
   nwfa=nwfa*1000000
   CALL mt19937_real2d(nifa(s:e,:))
   nifa=nifa*1000000
   CALL mt19937_real1d(nwfa2d(s:e))
   nwfa2d=61817
   CALL mt19937_real1d(nifa2d(s:e))
   nifa2d=0
   CALL mt19937_real2d(tgrs(s:e,:))
   tgrs=tgrs + 300.0
   CALL mt19937_real2d(prsl(s:e,:))
   prsl=prsl * 200000
   CALL mt19937_real2d(phil(s:e,:))
   phil=phil*800000
   CALL mt19937_real1d(area(s:e))
   area=2000000000.0
   CALL mt19937_real3d(diag3d(s:e,:,:))
   CALL mt19937_real2d(omega(s:e,:))
   omega=omega-1.0
   CALL mt19937_real1d(prcp(s:e))
   prcp=0
   CALL mt19937_real1d(rain(s:e))
   rain=0
   CALL mt19937_real1d(graupel(s:e))
   graupel=0
   CALL mt19937_real1d(ice(s:e))
   ice=0
   CALL mt19937_real1d(snow(s:e))
   snow=0
   CALL mt19937_real1d(sr(s:e))
   sr=0
   CALL mt19937_real2d(refl_10cm(s:e,:))
   refl_10cm=-35
   CALL mt19937_real2d(phii(s:e,:))
   phii=phii*800000
   CALL mt19937_real2d(spp_wts_mp(s:e,:))
   CALL mt19937_real1d(spp_prt_list(s:e))
   CALL mt19937_real1d(spp_stddev_cutoff(s:e))

   !MOHAMMAD
   !CALL mt19937_real1d(islmsk(s:e))
   islmsk=0
   CALL mt19937_real3d(aerfld(s:e,:,:))   
   CALL mt19937_real1d(max_hail_diam_sfc(s:e))
   max_hail_diam_sfc=0  
   CALL mt19937_real2d(pfi_lsan(s:e,:))
   pfi_lsan=0
   CALL mt19937_real2d(pfl_lsan(s:e,:))
   pfl_lsan=0 
         
   spp_var_list(:)=''
   !=============================================================

#ifdef MPI
   CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

   CALL SYSTEM_CLOCK (count_rate=count_rate)
   CALL SYSTEM_CLOCK (count=count_start)

   CALL SYSTEM_CLOCK (count=count_end)
   elapsed = REAL (count_end - count_start) / REAL (count_rate)
   PRINT*, "Finished copying data in =", elapsed  
   PRINT*

#ifdef MPI
   CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

   !--- Print state
   CALL print_state("Input state",   &
       ext_ndiag3d,      &
       spechum,          &
       qc,               &
       qr,               &
       qi,               &
       qs,               &
       qg,               &
       ni,               &
       nr,               &
       nc,               &
       nwfa,             &
       nifa,             &
       nwfa2d,           &
       nifa2d,           &
       tgrs,             &
       prsl,             &
       phil,             &
       area,             &
       diag3d,           &
       omega,            &
       prcp,             &
       rain,             &
       graupel,          &
       ice,              &
       snow,             &
       sr,               &
       refl_10cm,        &
       phii,             &
       spp_wts_mp,       &
       spp_prt_list,     &
       spp_stddev_cutoff &
   )
   !-------------

   PRINT*, "Calling init"
   CALL mp_thompson_init(ncol, nlev, con_g, con_rd, con_eps,     &
                        restart, imp_physics,                    &
                        imp_physics_thompson, convert_dry_rho,   &
                        spechum, qc, qr, qi, qs, qg, ni, nr,     &
                        is_aerosol_aware,  merra2_aerosol_aware, &
                        nc, nwfa2d, nifa2d,                      &
                        nwfa, nifa, tgrs, prsl, phil, area,      &
                        aerfld, mpicomm, mpirank, mpiroot,       &
                        threads, ext_diag, diag3d,               &
                        is_initialized, errmsg, errflg)                   

#ifdef MPI
   CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

   PRINT*, "Calling run"
   CALL SYSTEM_CLOCK (count_rate=count_rate)
   CALL SYSTEM_CLOCK (count=count_start)

   ncol_per_thread = (ncol / n_omp_threads)
   ncol_rem = n_omp_threads - mod(ncol, n_omp_threads)
!$omp parallel do private(tid,s,e)
   DO tid = 0, n_omp_threads - 1

       !--start and end--------
       if(tid <= ncol_rem) then
         s = tid * ncol_per_thread + 1
       else
         s = ncol_rem * ncol_per_thread + &
             (tid - ncol_rem) * (ncol_per_thread + 1) + 1
       endif
       if(tid < ncol_rem) then
         e = s + ncol_per_thread - 1
       else
         e = s + ncol_per_thread
       endif
       e = MIN(e, ncol)
       PRINT*, tid, s, ":", e, "=", e - s + 1
       !------------------------

       CALL mp_thompson_run(e-s+1, nlev, con_g, con_rd,        &
                            con_eps, convert_dry_rho,            &
                            spechum(s:e,:), qc(s:e,:), qr(s:e,:), qi(s:e,:), qs(s:e,:), qg(s:e,:), ni(s:e,:), nr(s:e,:), &
                            is_aerosol_aware, merra2_aerosol_aware, nc(s:e,:), nwfa(s:e,:), nifa(s:e,:),    &
                            nwfa2d(s:e), nifa2d(s:e), aero_ind_fdb,        &
                            tgrs(s:e,:), prsl(s:e,:), phii(s:e,:), omega(s:e,:),             &
                            sedi_semi, decfl, islmsk, dtp, dt_inner,     & 
                            first_time_step, istep, nsteps,      &
                            prcp(s:e), rain(s:e), graupel(s:e), ice(s:e), snow(s:e), sr(s:e),  &
                            refl_10cm(s:e,:), fullradar_diag, max_hail_diam_sfc, do_radar_ref, aerfld,  &
                            mpicomm, mpirank, mpiroot, blkno,    &
                            ext_diag, diag3d(s:e,:,:), reset_diag3d,      &
                            spp_wts_mp(s:e,:), spp_mp, n_var_spp,       &
                            spp_prt_list(s:e), spp_var_list(s:e),          &
                            spp_stddev_cutoff(s:e),                   &
                            cplchm, pfi_lsan, pfl_lsan,              &
                            is_initialized, errmsg, errflg)

   ENDDO
!$omp end parallel do

   CALL SYSTEM_CLOCK (count=count_end)
   elapsed = REAL (count_end - count_start) / REAL (count_rate)
   PRINT*
   PRINT*, "Finished executing kernel in =", elapsed  
   PRINT*

#ifdef MPI
   CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

   PRINT*, "Calling finalize"
   CALL mp_thompson_finalize(is_initialized, errmsg, errflg)

   !--- Print state
   CALL print_state("Output state",   &
       ext_ndiag3d,      &
       spechum,          &
       qc,               &
       qr,               &
       qi,               &
       qs,               &
       qg,               &
       ni,               &
       nr,               &
       nc,               &
       nwfa,             &
       nifa,             &
       nwfa2d,           &
       nifa2d,           &
       tgrs,             &
       prsl,             &
       phil,             &
       area,             &
       diag3d,           &
       omega,            &
       prcp,             &
       rain,             &
       graupel,          &
       ice,              &
       snow,             &
       sr,               &
       refl_10cm,        &
       phii,             &
       spp_wts_mp,       &
       spp_prt_list,     &
       spp_stddev_cutoff &
   )
   !-------------

#ifdef MPI
   CALL MPI_Finalize(ierror)
#endif

end program test_thompson

