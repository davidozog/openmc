module mic
  
  implicit none

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
