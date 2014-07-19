module mic
  
  implicit none

  type MICNuclide
    integer :: base_idx
    integer :: n_grid
    logical :: fissionable
    integer :: n_reaction
    integer :: index_fission(1)
    real(8) :: Q_value
    logical :: urr
  end type MICNuclide

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
  type(MICNuclide), allocatable, target :: mic_nuclides(:)
  integer, save ::  mic_n_nuclides_total 
  integer, save ::  mic_n_grid
  integer, save ::  mic_work
  public mic_materials, mic_n_nuclides, mic_grid_index, mic_energy, &
         mic_total, mic_elastic, mic_fission, mic_nu_fission, mic_absorption, &
         mic_heating, mic_nuclides, mic_n_nuclides_total, mic_n_grid, mic_work

!dir$ attributes offload:mic :: mic_materials, mic_n_nuclides, mic_grid_index, &
!dir$    mic_energy, mic_total, mic_elastic, mic_fission, mic_nu_fission,      &
!dir$    mic_absorption, mic_heating, mic_nuclides, mic_n_nuclides_total,      &
!dir$    mic_n_grid, mic_work
!dir$ attributes align:64 :: mic_materials, mic_n_nuclides, mic_grid_index,    &
!dir$    mic_energy, mic_total, mic_elastic, mic_fission, mic_nu_fission,      &
!dir$    mic_absorption, mic_heating, mic_nuclides, mic_n_nuclides_total,      &
!dir$    mic_n_grid, mic_work

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
