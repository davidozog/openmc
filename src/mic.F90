module mic
  
  implicit none

  integer, allocatable, save :: mic_materials(:)
  integer, allocatable, save :: mic_n_nuclides(:)
  integer, allocatable, save :: mic_base_nuclides(:)
  integer, allocatable, save :: mic_nuclides(:)
  public mic_materials, mic_n_nuclides, mic_base_nuclides, mic_nuclides

!dir$ attributes offload:mic :: mic_materials, mic_n_nuclides, &
!dir$                           mic_base_nuclides, mic_nuclides
!dir$ attributes align:64 :: mic_materials, mic_n_nuclides, mic_base_nuclides, &
!dir$                        mic_nuclides

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
