! These variables store particles for vector processing assistance.
! See Brown and Martin "MONTE CARLO METHODS FOR RADIATION TRANSPORT ANALYSIS
! ON VECTOR COMPUTERS" (1984) for details.

module mic
  
! TODO:
! use ace_header, only:NuclideMicroXS, MaterialMacroXS

  implicit none

! type MICNuclide
!   integer :: base_idx
!   integer :: n_grid
!   logical :: fissionable
!   integer :: n_reaction
!   integer :: index_fission(1)
!   real(8) :: Q_value
!   logical :: urr
! end type MICNuclide

!===============================================================================
! BANKEDPARTICLE constains the necessary state variables of a particle that
! needs to performance cross section lookups (for vectorization purposes)
!===============================================================================
!  type BankedParticle
!    ! Basic data
!    integer(8) :: id                  ! Unique ID
!    integer    :: type                ! Particle type (n, p, e, etc)
!    integer    :: material            ! index for current material
!    real(8)    :: E                   ! energy
!    integer    :: energy_index        ! unionized energy index
!    logical    :: check_sab           ! whether this material has S(a,b) tables
!    integer    :: n_nuclides          ! number of nuclides
!    integer    :: nuclides(MAX_NUCLIDES) ! Indices of the nuclides
!    real(8)    :: atom_density(MAX_NUCLIDES)  ! atom density of a nuclide
!  end type BankedParticle

  integer :: n_bank_xs = 0   ! Thread-private number of banked particles for 
                             ! cross section lookup

! type(BankedParticle), allocatable, target :: xs_bank(:)

! Banked Particle (bp) thread-private (tp) arrays
  integer(8), allocatable :: bp_tp_id(:)             ! Unique ID
  integer, allocatable    :: bp_tp_type(:)           ! Particle type (n, p, e, etc)
  integer, allocatable    :: bp_tp_material(:)       ! index for current material
  real(8), allocatable    :: bp_tp_E(:)              ! energy
  integer, allocatable    :: bp_tp_energy_index(:)   ! unionized energy index
  logical, allocatable    :: bp_tp_check_sab(:)      ! whether this material has S(a,b) tables
  integer, allocatable    :: bp_tp_n_nuclides(:)     ! number of nuclides
  integer, allocatable    :: bp_tp_nuclides(:,:)     ! Indices of the nuclides
  real(8), allocatable    :: bp_tp_atom_density(:,:) ! atom density of a nuclide

! Banked Particle (bp) arrays
  integer(8), allocatable :: bp_id(:)             ! Unique ID
  integer, allocatable    :: bp_type(:)           ! Particle type (n, p, e, etc)
  integer, allocatable    :: bp_material(:)       ! index for current material
  real(8), allocatable    :: bp_E(:)              ! energy
  integer, allocatable    :: bp_energy_index(:)   ! unionized energy index
  logical, allocatable    :: bp_check_sab(:)      ! whether this material has S(a,b) tables
  integer, allocatable    :: bp_n_nuclides(:)     ! number of nuclides
  integer, allocatable    :: bp_nuclides(:,:)     ! Indices of the nuclides
  real(8), allocatable    :: bp_atom_density(:,:) ! atom density of a nuclide
!$omp threadprivate(n_bank_xs, bp_tp_id, bp_tp_type, bp_tp_material, bp_tp_E, &
!$omp&              bp_tp_energy_index, bp_tp_check_sab, bp_tp_n_nuclides,    &
!$omp&              bp_tp_nuclides, bp_tp_atom_density)



  integer, allocatable, save :: mic_nuc_base_idx(:)
  integer, allocatable, save :: mic_nuc_Q_value(:)
! integer, allocatable :: mic_nuc_n_grid(:)
! logical, allocatable :: mic_nuc_fissionable(:)
! integer, allocatable :: mic_nuc_n_reaction(:)
! integer, allocatable :: mic_nuc_index_fission(:)
! real(8), allocatable :: mic_nuc_Q_value(:)
! logical, allocatable :: mic_nuc_urr(:)

  integer, allocatable, save :: mic_materials(:)
  integer, allocatable, save :: mic_n_nuclides(:)
  integer, allocatable, save :: mic_grid_index(:)
  real(8), allocatable, save :: mic_energy(:)
  real(8), allocatable, save :: mic_total(:)
  real(8), allocatable, save :: mic_elastic(:)
  real(8), allocatable, save :: mic_fission(:)
  real(8), allocatable, save :: mic_nu_fission(:)
  real(8), allocatable, save :: mic_absorption(:)
  real(8), allocatable, save :: mic_heating(:)
! type(MICNuclide), allocatable, target :: mic_nuclides(:)

!   type(BankedParticle), intent(inout) :: master_xs_bank(:)
!!dir$ attributes align : 64 :: mymicro_xs, mymaterial_xs

! type(NuclideMicroXS) :: mymicro_xs(n_nuclides_total, work)  ! Cache for each nuclide
  integer,allocatable :: mic_micro_index_grid(:,:)         ! index on nuclide energy grid
  integer,allocatable :: mic_micro_index_temp(:,:)         ! temperature index for nuclide
  real(8),allocatable :: mic_micro_last_E(:,:)             ! last evaluated energy
  real(8),allocatable :: mic_micro_interp_factor(:,:)      ! interpolation factor on nuc. energy grid
  real(8),allocatable :: mic_micro_total(:,:)              ! microscropic total xs
  real(8),allocatable :: mic_micro_elastic(:,:)            ! microscopic elastic scattering xs
  real(8),allocatable :: mic_micro_absorption(:,:)         ! microscopic absorption xs
  real(8),allocatable :: mic_micro_fission(:,:)            ! microscopic fission xs
  real(8),allocatable :: mic_micro_nu_fission(:,:)         ! microscopic production xs
  real(8),allocatable :: mic_micro_kappa_fission(:,:)      ! microscopic energy-released from fission
  integer,allocatable :: mic_micro_index_sab(:,:)          !(work) index in sab_tables (zero means no table)
  integer,allocatable :: mic_micro_last_index_sab(:,:)     ! index in sab_tables last used by this nuclide
  real(8),allocatable :: mic_micro_elastic_sab(:,:)        ! microscopic elastic scattering on S(a,b) table
  logical,allocatable :: mic_micro_use_ptable(:,:)         ! in URR range with probability tables?

! type(MaterialMacroXS) :: mymaterial_xs(work)  ! Cache for each nuclide
  real(8), allocatable :: mic_mat_total(:)
  real(8), allocatable :: mic_mat_elastic(:)        ! macroscopic elastic scattering xs
  real(8), allocatable :: mic_mat_absorption(:)     ! macroscopic absorption xs
  real(8), allocatable :: mic_mat_fission(:)        ! macroscopic fission xs
  real(8), allocatable :: mic_mat_nu_fission(:)     ! macroscopic production xs
  real(8), allocatable :: mic_mat_kappa_fission(:)  ! macroscopic energy-released from fission

!  TODO:
!  type(NuclideMicroXS), allocatable :: micro_xs(:)  ! Cache for each nuclide
!  type(MaterialMacroXS)             :: material_xs  ! Cache for current material

  integer, save ::  mic_n_nuclides_total 
  integer, save ::  mic_n_grid
  integer, save ::  mic_work
  public mic_materials, mic_n_nuclides, mic_grid_index, mic_energy, &
         mic_total, mic_elastic, mic_fission, mic_nu_fission, mic_absorption, &
         mic_heating, mic_nuc_base_idx, mic_n_nuclides_total, mic_n_grid,     &
         mic_nuc_Q_value, mic_work, bp_tp_id, bp_tp_type, bp_tp_material,     &
         bp_tp_E, bp_tp_energy_index, bp_tp_check_sab, bp_tp_n_nuclides,      &
         bp_tp_nuclides, bp_tp_atom_density

!dir$ attributes offload:mic :: mic_materials, mic_n_nuclides, mic_grid_index, &
!dir$    mic_energy, mic_total, mic_elastic, mic_fission, mic_nu_fission,      &
!dir$    mic_absorption, mic_heating, mic_nuc_base_idx, mic_n_nuclides_total,      &
!dir$    mic_n_grid, mic_work, mic_micro_index_grid, mic_micro_index_temp,     &
!dir$    mic_micro_last_E, mic_micro_interp_factor, mic_micro_total,           &
!dir$    mic_micro_elastic, mic_micro_absorption, mic_micro_fission,           &
!dir$    mic_micro_nu_fission, mic_micro_kappa_fission, mic_micro_index_sab,   &
!dir$    mic_micro_last_index_sab, mic_micro_elastic_sab, mic_micro_use_ptable,&
!dir$    mic_mat_total, mic_mat_elastic, mic_mat_absorption, mic_mat_fission,  &
!dir$    mic_mat_nu_fission, mic_mat_kappa_fission, mic_nuc_Q_value,           &
!dir$    bp_id, bp_type, bp_material, bp_E, bp_energy_index,                   &
!dir$    bp_check_sab, bp_n_nuclides, bp_nuclides, bp_atom_density,            &
!dir$    bp_tp_id, bp_tp_type, bp_tp_material, bp_tp_E, bp_tp_energy_index,    &
!dir$    bp_tp_check_sab, bp_tp_n_nuclides, bp_tp_nuclides, bp_tp_atom_density

!dir$ attributes align:64 :: mic_materials, mic_n_nuclides, mic_grid_index,    &
!dir$    mic_energy, mic_total, mic_elastic, mic_fission, mic_nu_fission,      &
!dir$    mic_absorption, mic_heating, mic_nuc_base_idx, mic_n_nuclides_total,      &
!dir$    mic_n_grid, mic_work, mic_micro_index_grid, mic_micro_index_temp,     &
!dir$    mic_micro_last_E, mic_micro_interp_factor, mic_micro_total,           &
!dir$    mic_micro_elastic, mic_micro_absorption, mic_micro_fission,           &
!dir$    mic_micro_nu_fission, mic_micro_kappa_fission, mic_micro_index_sab,   &
!dir$    mic_micro_last_index_sab, mic_micro_elastic_sab, mic_micro_use_ptable,&
!dir$    mic_mat_total, mic_mat_elastic, mic_mat_absorption, mic_mat_fission,  &
!dir$    mic_mat_nu_fission, mic_mat_kappa_fission, mic_nuc_Q_value,           &
!dir$    bp_id, bp_type, bp_material, bp_E, bp_energy_index,                   &
!dir$    bp_check_sab, bp_n_nuclides, bp_nuclides, bp_atom_density,            &
!dir$    bp_tp_id, bp_tp_type, bp_tp_material, bp_tp_E, bp_tp_energy_index,    &
!dir$    bp_tp_check_sab, bp_tp_n_nuclides, bp_tp_nuclides, bp_tp_atom_density

contains 

  subroutine mic_func()

    integer :: total_xs      ! number of particles in bank
    integer :: pp            ! particle index

    total_xs = 100
!dir$ offload begin target(mic:0)
!$omp parallel do
    do pp = 1, total_xs
      print *, pp
    end do
!dir$ end offload 
    
  end subroutine mic_func

end module mic
