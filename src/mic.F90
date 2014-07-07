module mic
  
  implicit none

  type MICNuclide
    integer :: base_idx
    integer :: n_grid
  end type MICNuclide

  integer, allocatable, save :: mic_materials(:)
  integer, allocatable, save :: mic_n_nuclides(:)
  integer, allocatable, save :: mic_grid_index(:)
  integer, allocatable, save :: mic_energy(:)
  integer, allocatable, save :: mic_total(:)
  integer, allocatable, save :: mic_elastic(:)
  integer, allocatable, save :: mic_fission(:)
  integer, allocatable, save :: mic_nu_fission(:)
  integer, allocatable, save :: mic_absorption(:)
  integer, allocatable, save :: mic_heating(:)
  type(MICNuclide), allocatable, target :: mic_nuclides(:)
  integer, save ::  mic_n_nuclides_total 
  public mic_materials, mic_n_nuclides, mic_grid_index, mic_energy, &
         mic_total, mic_elastic, mic_fission, mic_nu_fission, mic_absorption, &
         mic_heating, mic_nuclides, mic_n_nuclides_total

!dir$ attributes offload:mic :: mic_materials, mic_n_nuclides, mic_grid_index, &
!dir$    mic_energy, mic_total, mic_elastic, mic_fission, mic_nu_fission,      &
!dir$    mic_absorption, mic_heating, mic_nuclides, mic_n_nuclides_total
!dir$ attributes align:64 :: mic_materials, mic_n_nuclides, mic_grid_index,    &
!dir$    mic_energy, mic_total, mic_elastic, mic_fission, mic_nu_fission,      &
!dir$    mic_absorption, mic_heating, mic_nuclides, mic_n_nuclides_total 

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
