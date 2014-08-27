module cross_section

  use ace_header,      only: Nuclide, SAlphaBeta, Reaction, UrrData
  use constants
  use error,           only: fatal_error
  use fission,         only: nu_total
  use global
  use material_header, only: Material
  use particle_header, only: Particle
  use random_lcg,      only: prn
  use search,          only: binary_search
  use mic
  use omp_lib

  implicit none
  save
  
  integer :: union_grid_index

!!dir$ attributes align:64 :: master_xs_bank
!  type(bankedparticle), allocatable, target :: master_xs_bank(:)
!  public :: master_xs_bank

!$omp threadprivate(union_grid_index)

contains

!===============================================================================
! XS_BANKS places particle in a bank for vectorized cross section lookups later
!===============================================================================

  subroutine bank_xs(p)

    integer :: current_bank_slot
    integer :: i
    type(Particle), intent(inout) :: p
    type(Material), pointer :: mat => null() ! current material
    
    current_bank_slot = n_bank_xs + 1
    bp_tp_id(current_bank_slot) = p % id
!   bp_tp_type(current_bank_slot) = p % type
!   bp_tp_material(current_bank_slot) = p % material
    bp_tp_E(current_bank_slot) = p % E

    mat => materials(p % material)
!   bp_tp_check_sab(current_bank_slot) = (mat % n_sab > 0)
    bp_tp_n_nuclides(current_bank_slot) = mat % n_nuclides

    if (mat % n_nuclides .ge. MAX_NUCLIDES) print *, "MAXED OUT NUCS!"

    do i=1, mat % n_nuclides
      bp_tp_nuclides(i, current_bank_slot) = mat % nuclide(i)
      bp_tp_atom_density(i, current_bank_slot) = mat % atom_density(i)
    end do

    ! increment number of bank sites
#ifdef _OPENMP
    n_bank_xs = min(n_bank_xs + 1, work/n_threads+2)
    if (n_bank_xs .eq. work/n_threads+2) print *, "MAXED OUT!", thread_id
#else
    n_bank_xs = min(int(n_bank_xs + 1, 8), work+1)
#endif

    if (p % material == MATERIAL_VOID) then
      bp_tp_energy_index(current_bank_slot) = -1
      print *, "WARNING: NOT SURE IF THIS IS COOL..."
    else if (grid_method == GRID_UNION) then
      bp_tp_energy_index(current_bank_slot) = &
              find_energy_index(bp_tp_E(current_bank_slot))
    else
      message = "NON UNION MIC NOT IMPLEMENTED YET."
!!!     if (E < nuc % energy(1)) then
!!!       i_grid = 1
!!!     elseif (E > nuc % energy(nuc % n_grid)) then
!!!       i_grid = nuc % n_grid - 1
!!!     else
!!!       i_grid = binary_search(nuc % energy, nuc % n_grid, E)
!!!     end if
      call fatal_error()
    end if

  end subroutine bank_xs

!===============================================================================
! CALCULATE_BANK_XS does the same thing as calculate_xs, but vectorized for MIC
!===============================================================================

  subroutine calculate_bank_xs()
    integer :: i             ! loop index over nuclides
    integer :: p_id          ! particle ID 
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in mat % i_sab_nuclides
    integer :: pp            ! particle index
    integer :: total_xs      ! number of particles in bank
    real(8) :: atom_density  ! atom density of a nuclide
    logical :: check_sab     ! should we check for S(a,b) table?
    type(Material), pointer, save :: mat => null() ! current material

    integer :: E_idx         ! index into grid_index
    integer :: i_grid        ! index on nuclide energy grid
    real(8) :: E             ! energy
    real(8) :: f             ! interp factor on nuclide energy grid
    real(8) :: wtime

    print *, "calculate_bank_xs() ..."

    ! Set all material macroscopic cross sections to zero
    do i = 1, work
      mic_mat_total(i)      = ZERO
      mic_mat_elastic(i)    = ZERO
      mic_mat_absorption(i) = ZERO
      mic_mat_fission(i)    = ZERO
      mic_mat_nu_fission(i) = ZERO
      mic_mat_kappa_fission(i)  = ZERO
    end do

    total_xs = n_bank_xs
!    print *, "n_grid:", n_grid
    print *, "total_xs is ", total_xs
!    print *, "bank:", xs_bank

!    do i = 1, work
!      do j = 1, n_nuclides_total
!        mymicro_xs(j, i) = micro_xs(j)
!!     print *, "MYMICRO-host:", mymicro_xs(i)
!      end do
!    end do

!   print *, "DOING:", total_xs

    wtime = omp_get_wtime()

!$omp parallel do private(atom_density,i_nuclide,i_sab, p_id, i_grid, E, E_idx, f, i) &
!$omp&            schedule(guided)
    do pp = 1, total_xs
!     print *, "hi:", omp_get_thread_num(), master_xs_bank(pp)

!     write(*,*)'Thread id: ', OMP_GET_THREAD_NUM()

!     j = 1
      p_id = bp_id(pp)

!   do i=1, mic_n_nuclides_total
!     print *, "MYMIRCRO-mic:", mymicro_xs(i)
!   end do

!!dir$ noprefetch mic_micro_index_grid
!!dir$ prefetch bp_n_nuclides, mic_grid_index, mic_nuc_base_idx
!!dir$ prefetch bp_nuclides, bp_energy_index, bp_E

!NOTE: might need a "simd" pragma here
!dir$ ivdep
      NUCLIDES_IN_MATERIAL_LOOP: do i = 1, bp_n_nuclides(pp)
        ! ======================================================================
        ! CHECK FOR S(A,B) TABLE
        i_sab = 0
        ! Check if this nuclide matches one of the S(a,b) tables specified --
        ! this relies on i_sab_nuclides being in sorted order
!        if (master_xs_bank(pp) % check_sab) then
!          if (i == mat % i_sab_nuclides(j)) then
!            ! Get index in sab_tables
!            i_sab = mat % i_sab_tables(j)
!
!            ! If particle energy is greater than highest energy for the
!            ! S(a,v) table, don't use the S(a,b) table
!            if (master_xs_bank(pp) % E > sab_tables(i_sab) % threshold_inelastic) i_sab = 0
!
!            ! Increment position in i_sab_nuclides
!            j = j + 1
!
!            ! Don't check for S(a,b) tables if there are no more left
!            if (j > mat % n_sab) check_sab = .false.
!          end if
!        end if

!       TODO:  Do I need this?  Why MAX used?
        i_nuclide = bp_nuclides(i, pp)
        ! Calculate microscopic cross section for this nuclide
!       if (master_xs_bank(pp) % E /= mymicro_xs(i_nuclide, p_id) % last_E) then
!         call calculate_nuclide_xs(i_nuclide, i_sab, master_xs_bank(pp) % E, &
!             mymicro_xs, master_xs_bank(pp) % energy_index, p_id)
!       else if (i_sab /= mymicro_xs(i_nuclide, p_id) % last_index_sab) then
!           call calculate_nuclide_xs(i_nuclide, i_sab, master_xs_bank(pp) % E, &
!               mymicro_xs, master_xs_bank(pp) % energy_index, p_id)
!       end if

!   nuc => mic_nuclides(i_nuclide)
    i_grid = bp_energy_index(pp)
    E = bp_E(pp)
    !E_idx = mic_grid_index((i_nuclide-1)*mic_n_grid + i_grid)
    E_idx = mic_nuc_base_idx(i_nuclide) + mic_grid_index((i_nuclide-1)*mic_n_grid + i_grid)

    ! check for rare case where two energy points are the same
    !if (mic_energy(E_idx) == mic_energy(E_idx+1)) E_idx = E_idx + 1
    f = (E - mic_energy(E_idx))/(mic_energy(E_idx+1) - mic_energy(E_idx))

    mic_micro_index_grid(i_nuclide, p_id)     = E_idx
    mic_micro_interp_factor(i_nuclide, p_id) = f

    ! Initialize nuclide cross-sections to zero
    mic_micro_fission(i_nuclide, p_id) = ZERO
    mic_micro_nu_fission(i_nuclide, p_id) = ZERO
    mic_micro_kappa_fission(i_nuclide, p_id) = ZERO

    ! Calculate microscopic nuclide total cross section
    mic_micro_total(i_nuclide, p_id) = (ONE - f) * mic_total(E_idx) &
          + f * mic_total(E_idx+1)

    ! Calculate microscopic nuclide total cross section
    mic_micro_elastic(i_nuclide, p_id) = (ONE - f) * mic_elastic(E_idx) &
          + f * mic_elastic(E_idx+1)

    ! Calculate microscopic nuclide absorption cross section
    mic_micro_absorption(i_nuclide, p_id) = (ONE - f) * mic_absorption( &
          E_idx) + f * mic_absorption(E_idx+1)

!   TODO: Could feasibly remove conditionals!
!   if (nuc % fissionable) then
      ! Calculate microscopic nuclide total cross section
      mic_micro_fission(i_nuclide, p_id) = (ONE - f) * mic_fission(E_idx) &
           + f * mic_fission(E_idx+1)

      ! Calculate microscopic nuclide nu-fission cross section
      mic_micro_nu_fission(i_nuclide, p_id) = (ONE - f) * mic_nu_fission( &
           E_idx) + f * mic_nu_fission(E_idx+1)
          
      ! Calculate microscopic nuclide kappa-fission cross section
      ! The ENDF standard (ENDF-102) states that MT 18 stores
      ! the fission energy as the Q_value (fission(1))
      mic_micro_kappa_fission(i_nuclide, p_id) = &
           mic_nuc_Q_value(i_nuclide) * mic_micro_fission(i_nuclide, p_id)
!   end if

    ! If there is S(a,b) data for this nuclide, we need to do a few
    ! things. Since the total cross section was based on non-S(a,b) data, we
    ! need to correct it by subtracting the non-S(a,b) elastic cross section and
    ! then add back in the calculated S(a,b) elastic+inelastic cross section.

!    if (i_sab > 0) then
!!     print *, "SAB NOT IMPLEMENTED YET."
!      stop 1
!      !call calculate_sab_xs(i_nuclide, i_sab, E)
!    end if

    ! if the particle is in the unresolved resonance range and there are
    ! probability tables, we need to determine cross sections from the table
 
    !NOTE: URR XS NOT SUPPORTED YET ON MIC
    !   if (urr_ptables_on .and. nuc % urr_present) then
    !   if (urr_ptables_on .and. nuc % urr_present) then
    !     if (E > nuc % urr_data % energy(1) .and. &
    !          E < nuc % urr_data % energy(nuc % urr_data % n_energy)) then
    !       call calculate_urr_xs(i_nuclide, E)
    !     end if
    !   end if

    mic_micro_last_E(i_nuclide, p_id) = E
    mic_micro_last_index_sab(i_nuclide, p_id) = i_sab

        ! Copy atom density of nuclide in material
        atom_density = bp_atom_density(i, pp)

        ! Add contributions to material macroscopic total cross section
        mic_mat_total(p_id) = mic_mat_total(p_id) + &
             atom_density * mic_micro_total(i_nuclide, p_id) 

        ! Add contributions to material macroscopic scattering cross section
        mic_mat_elastic(p_id) = mic_mat_elastic(p_id) + &
             atom_density * mic_micro_elastic(i_nuclide, p_id)

        ! Add contributions to material macroscopic absorption cross section
        mic_mat_absorption(p_id) = mic_mat_absorption(p_id) + & 
             atom_density * mic_micro_absorption(i_nuclide, p_id) 

        ! Add contributions to material macroscopic fission cross section
        mic_mat_fission(p_id) = mic_mat_fission(p_id) + &
             atom_density * mic_micro_fission(i_nuclide, p_id)

        ! Add contributions to material macroscopic nu-fission cross section
        mic_mat_nu_fission(p_id) = mic_mat_nu_fission(p_id) + &
             atom_density * mic_micro_nu_fission(i_nuclide, p_id)
             
        ! Add contributions to material macroscopic energy release from fission
        mic_mat_kappa_fission(p_id) = mic_mat_kappa_fission(p_id) + &
             atom_density * mic_micro_kappa_fission(i_nuclide, p_id) 

!      print *, "MYMIRCRO-after-call:", mymicro_xs(i_nuclide, p_id) % total, &
!      "i_nuclide:", i_nuclide

!       print *, "particle:", master_xs_bank(pp) % id, "nuc", i, &
!                "material_xs:", mymaterial_xs(p_id)
!      print *, "micro_xs:", mymicro_xs(i_nuclide, p_id)
      
      end do NUCLIDES_IN_MATERIAL_LOOP
!     print *, "pp=", pp, "i=", i

    end do 
!$omp end parallel do
    wtime = omp_get_wtime() - wtime

    print *, "time:", wtime

!  print *, "mymaterials for each particle ID:"
!  do i = 1, total_xs
!    write(*,*), i, mic_mat_total(i), mic_mat_elastic(i), mic_mat_absorption(i),   &
!         mic_mat_fission(i), mic_mat_nu_fission(i)!, mic_mat_kappa_fission(i)
!  end do
!
!    do pp = 1, total_xs
!
!      print *, "here" 
!
!      ! Exit subroutine if material is void
!      !if (xs_bank(pp) % material == MATERIAL_VOID) return
!
!      mat => materials(xs_bank(pp) % material)
!
!      ! Find energy index on unionized grid
!      !if (grid_method == GRID_UNION) call find_energy_index(xs_bank(pp) % E)
!
!      ! Determine if this material has S(a,b) tables
!      check_sab = (mat % n_sab > 0)
!
!      ! Initialize position in i_sab_nuclides
!      j = 1
!
!      ! Add contribution from each nuclide in material
!      do i = 1, mat % n_nuclides
!        ! ======================================================================
!        ! CHECK FOR S(A,B) TABLE
!
!!       i_sab = 0
!
!        ! Check if this nuclide matches one of the S(a,b) tables specified --
!        ! this relies on i_sab_nuclides being in sorted order
!        if (check_sab) then
!          if (i == mat % i_sab_nuclides(j)) then
!            ! Get index in sab_tables
!            i_sab = mat % i_sab_tables(j)
!
!            ! If particle energy is greater than highest energy for the
!            ! S(a,v) table, don't use the S(a,b) table
!            if (xs_bank(pp) % E > sab_tables(i_sab) % threshold_inelastic) i_sab = 0
!
!            ! Increment position in i_sab_nuclides
!            j = j + 1
!
!            ! Don't check for S(a,b) tables if there are no more left
!            if (j > mat % n_sab) check_sab = .false.
!          end if
!        end if
!
!        ! ======================================================================
!        ! CALCULATE MICROSCOPIC CROSS SECTION
!
!        ! Determine microscopic cross sections for this nuclide
!!       i_nuclide = mat % nuclide(i)
!
!        ! Calculate microscopic cross section for this nuclide
!        !if (xs_bank(pp) % E /= micro_xs(i_nuclide) % last_E) then
!        !  call calculate_nuclide_xs(i_nuclide, i_sab, xs_bank(pp) % E)
!        !else if (i_sab /= micro_xs(i_nuclide) % last_index_sab) then
!        !  call calculate_nuclide_xs(i_nuclide, i_sab, xs_bank(pp) % E)
!        !end if
!
!        ! ======================================================================
!        ! ADD TO MACROSCOPIC CROSS SECTION
!
!        ! Copy atom density of nuclide in material
!!       atom_density = mat % atom_density(i)
!
!!       ! Add contributions to material macroscopic total cross section
!!       material_xs % total = material_xs % total + &
!!            atom_density * micro_xs(i_nuclide) % total
!
!!       ! Add contributions to material macroscopic scattering cross section
!!       material_xs % elastic = material_xs % elastic + &
!!            atom_density * micro_xs(i_nuclide) % elastic
!
!!       ! Add contributions to material macroscopic absorption cross section
!!       material_xs % absorption = material_xs % absorption + & 
!!            atom_density * micro_xs(i_nuclide) % absorption
!
!!       ! Add contributions to material macroscopic fission cross section
!!       material_xs % fission = material_xs % fission + &
!!            atom_density * micro_xs(i_nuclide) % fission
!
!!       ! Add contributions to material macroscopic nu-fission cross section
!!       material_xs % nu_fission = material_xs % nu_fission + &
!!            atom_density * micro_xs(i_nuclide) % nu_fission
!!            
!!       ! Add contributions to material macroscopic energy release from fission
!!       material_xs % kappa_fission = material_xs % kappa_fission + &
!!            atom_density * micro_xs(i_nuclide) % kappa_fission
!      end do
!    end do

!$omp parallel 
    n_bank_xs = 0
!$omp end parallel 

  end subroutine calculate_bank_xs

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle), intent(inout) :: p

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in mat % i_sab_nuclides
    real(8) :: atom_density  ! atom density of a nuclide
    logical :: check_sab     ! should we check for S(a,b) table?
    type(Material), pointer, save :: mat => null() ! current material
!$omp threadprivate(mat)

    ! Set all material macroscopic cross sections to zero
    material_xs % total      = ZERO
    material_xs % elastic    = ZERO
    material_xs % absorption = ZERO
    material_xs % fission    = ZERO
    material_xs % nu_fission = ZERO
    material_xs % kappa_fission  = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    mat => materials(p % material)

    ! Find energy index on unionized grid
    if (grid_method == GRID_UNION) union_grid_index = find_energy_index(p % E)

    ! Determine if this material has S(a,b) tables
    check_sab = (mat % n_sab > 0)

    ! Initialize position in i_sab_nuclides
    j = 1

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
      ! ========================================================================
      ! CHECK FOR S(A,B) TABLE

      i_sab = 0

      ! Check if this nuclide matches one of the S(a,b) tables specified -- this
      ! relies on i_sab_nuclides being in sorted order
!!      if (check_sab) then
!!        if (i == mat % i_sab_nuclides(j)) then
!!          ! Get index in sab_tables
!!          i_sab = mat % i_sab_tables(j)
!!
!!          ! If particle energy is greater than the highest energy for the S(a,b)
!!          ! table, don't use the S(a,b) table
!!          if (p % E > sab_tables(i_sab) % threshold_inelastic) i_sab = 0
!!
!!          ! Increment position in i_sab_nuclides
!!          j = j + 1
!!
!!          ! Don't check for S(a,b) tables if there are no more left
!!          if (j > mat % n_sab) check_sab = .false.
!!        end if
!!      end if

      ! ========================================================================
      ! CALCULATE MICROSCOPIC CROSS SECTION

      ! Determine microscopic cross sections for this nuclide
      i_nuclide = mat % nuclide(i)

      ! Calculate microscopic cross section for this nuclide
      if (p % E /= micro_xs(i_nuclide) % last_E) then
        call calculate_nuclide_xs(i_nuclide, i_sab, p % E)
      else if (i_sab /= micro_xs(i_nuclide) % last_index_sab) then
        call calculate_nuclide_xs(i_nuclide, i_sab, p % E)
      end if

      ! ========================================================================
      ! ADD TO MACROSCOPIC CROSS SECTION

      ! Copy atom density of nuclide in material
      atom_density = mat % atom_density(i)

      ! Add contributions to material macroscopic total cross section
      material_xs % total = material_xs % total + &
           atom_density * micro_xs(i_nuclide) % total

      ! Add contributions to material macroscopic scattering cross section
      material_xs % elastic = material_xs % elastic + &
           atom_density * micro_xs(i_nuclide) % elastic

      ! Add contributions to material macroscopic absorption cross section
      material_xs % absorption = material_xs % absorption + & 
           atom_density * micro_xs(i_nuclide) % absorption

      ! Add contributions to material macroscopic fission cross section
      material_xs % fission = material_xs % fission + &
           atom_density * micro_xs(i_nuclide) % fission

      ! Add contributions to material macroscopic nu-fission cross section
      material_xs % nu_fission = material_xs % nu_fission + &
           atom_density * micro_xs(i_nuclide) % nu_fission
           
      ! Add contributions to material macroscopic energy release from fission
      material_xs % kappa_fission = material_xs % kappa_fission + &
           atom_density * micro_xs(i_nuclide) % kappa_fission
    end do

  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(i_nuclide, i_sab, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy

    integer :: i_grid ! index on nuclide energy grid
    real(8) :: f      ! interp factor on nuclide energy grid
    type(Nuclide), pointer, save :: nuc => null()
!$omp threadprivate(nuc)

    ! Set pointer to nuclide
    nuc => nuclides(i_nuclide)

    ! Determine index on nuclide energy grid
    select case (grid_method)
    case (GRID_UNION)
      ! If we're using the unionized grid with pointers, finding the index on
      ! the nuclide energy grid is as simple as looking up the pointer

      i_grid = nuc % grid_index(union_grid_index)

    case (GRID_NUCLIDE)
      ! If we're not using the unionized grid, we have to do a binary search on
      ! the nuclide energy grid in order to determine which points to
      ! interpolate between

      if (E < nuc % energy(1)) then
        i_grid = 1
      elseif (E > nuc % energy(nuc % n_grid)) then
        i_grid = nuc % n_grid - 1
      else
        i_grid = binary_search(nuc % energy, nuc % n_grid, E)
      end if

    end select

    ! check for rare case where two energy points are the same
    if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1

    ! calculate interpolation factor
    f = (E - nuc%energy(i_grid))/(nuc%energy(i_grid+1) - nuc%energy(i_grid))

    micro_xs(i_nuclide) % index_grid    = i_grid
    micro_xs(i_nuclide) % interp_factor = f

    ! Initialize sab treatment to false
    micro_xs(i_nuclide) % index_sab   = NONE
    micro_xs(i_nuclide) % elastic_sab = ZERO
    micro_xs(i_nuclide) % use_ptable  = .false.

    ! Initialize nuclide cross-sections to zero
    micro_xs(i_nuclide) % fission    = ZERO
    micro_xs(i_nuclide) % nu_fission = ZERO
    micro_xs(i_nuclide) % kappa_fission  = ZERO

    ! Calculate microscopic nuclide total cross section
!!    micro_xs(i_nuclide) % total = (ONE - f) * nuc % total(i_grid) &
!!         + f * nuc % total(i_grid+1)
!!  TODO: FIX THIS WEIRDNESS!
    micro_xs(i_nuclide) % total = (ONE - f) * nuc % fission(i_grid) &
         + f * nuc % fission(i_grid+1)

    ! Calculate microscopic nuclide total cross section
    micro_xs(i_nuclide) % elastic = (ONE - f) * nuc % elastic(i_grid) &
         + f * nuc % elastic(i_grid+1)

    ! Calculate microscopic nuclide absorption cross section
    micro_xs(i_nuclide) % absorption = (ONE - f) * nuc % absorption( &
         i_grid) + f * nuc % absorption(i_grid+1)

!!    if (nuc % fissionable) then
!!      ! Calculate microscopic nuclide total cross section
!!      micro_xs(i_nuclide) % fission = (ONE - f) * nuc % fission(i_grid) &
!!           + f * nuc % fission(i_grid+1)
!!
!!      ! Calculate microscopic nuclide nu-fission cross section
!!      micro_xs(i_nuclide) % nu_fission = (ONE - f) * nuc % nu_fission( &
!!           i_grid) + f * nuc % nu_fission(i_grid+1)
!!           
!!      ! Calculate microscopic nuclide kappa-fission cross section
!!      ! The ENDF standard (ENDF-102) states that MT 18 stores
!!      ! the fission energy as the Q_value (fission(1))
!!      micro_xs(i_nuclide) % kappa_fission = &
!!           nuc % reactions(nuc % index_fission(1)) % Q_value * &
!!           micro_xs(i_nuclide) % fission
!!    end if

    ! If there is S(a,b) data for this nuclide, we need to do a few
    ! things. Since the total cross section was based on non-S(a,b) data, we
    ! need to correct it by subtracting the non-S(a,b) elastic cross section and
    ! then add back in the calculated S(a,b) elastic+inelastic cross section.

!!    if (i_sab > 0) call calculate_sab_xs(i_nuclide, i_sab, E)
!! 
!!   ! if the particle is in the unresolved resonance range and there are
!!   ! probability tables, we need to determine cross sections from the table
!! 
!!     if (urr_ptables_on .and. nuc % urr_present) then
!!       if (E > nuc % urr_data % energy(1) .and. &
!!            E < nuc % urr_data % energy(nuc % urr_data % n_energy)) then
!!         call calculate_urr_xs(i_nuclide, E)
!!       end if
!!     end if
!! 
     micro_xs(i_nuclide) % last_E = E
     micro_xs(i_nuclide) % last_index_sab = i_sab

  end subroutine calculate_nuclide_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace
! whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(i_nuclide, i_sab, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy

    integer :: i_grid    ! index on S(a,b) energy grid
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section
    type(SAlphaBeta), pointer, save :: sab => null()
!$omp threadprivate(sab)

    ! Set flag that S(a,b) treatment should be used for scattering
    micro_xs(i_nuclide) % index_sab = i_sab

    ! Get pointer to S(a,b) table
    sab => sab_tables(i_sab)

    ! Get index and interpolation factor for inelastic grid
    if (E < sab % inelastic_e_in(1)) then
      i_grid = 1
      f = ZERO
    else
      i_grid = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, E)
      f = (E - sab%inelastic_e_in(i_grid)) / & 
           (sab%inelastic_e_in(i_grid+1) - sab%inelastic_e_in(i_grid))
    end if

    ! Calculate S(a,b) inelastic scattering cross section
    inelastic = (ONE - f) * sab % inelastic_sigma(i_grid) + &
         f * sab % inelastic_sigma(i_grid + 1)

    ! Check for elastic data
    if (E < sab % threshold_elastic) then
      ! Determine whether elastic scattering is given in the coherent or
      ! incoherent approximation. For coherent, the cross section is
      ! represented as P/E whereas for incoherent, it is simply P

      if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
        if (E < sab % elastic_e_in(1)) then
          ! If energy is below that of the lowest Bragg peak, the elastic
          ! cross section will be zero
          elastic = ZERO
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, E)
          elastic = sab % elastic_P(i_grid) / E
        end if
      else
        ! Determine index on elastic energy grid
        if (E < sab % elastic_e_in(1)) then
          i_grid = 1
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, E)
        end if

        ! Get interpolation factor for elastic grid
        f = (E - sab%elastic_e_in(i_grid))/(sab%elastic_e_in(i_grid+1) - &
             sab%elastic_e_in(i_grid))

        ! Calculate S(a,b) elastic scattering cross section
        elastic = (ONE - f) * sab % elastic_P(i_grid) + &
             f * sab % elastic_P(i_grid + 1)
      end if
    else
      ! No elastic data
      elastic = ZERO
    end if

    ! Correct total and elastic cross sections
    micro_xs(i_nuclide) % total = micro_xs(i_nuclide) % total - &
         micro_xs(i_nuclide) % elastic + inelastic + elastic
    micro_xs(i_nuclide) % elastic = inelastic + elastic

    ! Store S(a,b) elastic cross section for sampling later
    micro_xs(i_nuclide) % elastic_sab = elastic

  end subroutine calculate_sab_xs

!===============================================================================
! CALCULATE_URR_XS determines cross sections in the unresolved resonance range
! from probability tables
!===============================================================================

  subroutine calculate_urr_xs(i_nuclide, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    real(8), intent(in) :: E         ! energy

    integer :: i_energy   ! index for energy
    integer :: i_table    ! index for table
    real(8) :: f          ! interpolation factor
    real(8) :: r          ! pseudo-random number
    real(8) :: elastic    ! elastic cross section
    real(8) :: capture    ! (n,gamma) cross section
    real(8) :: fission    ! fission cross section
    real(8) :: inelastic  ! inelastic cross section
    type(UrrData),  pointer, save :: urr => null()
    type(Nuclide),  pointer, save :: nuc => null()
    type(Reaction), pointer, save :: rxn => null()
!$omp threadprivate(urr, nuc, rxn)

    micro_xs(i_nuclide) % use_ptable = .true.

    ! get pointer to probability table
    nuc => nuclides(i_nuclide)
    urr => nuc % urr_data

    ! determine energy table
    i_energy = 1
    do
      if (E < urr % energy(i_energy + 1)) exit
      i_energy = i_energy + 1
    end do

    ! determine interpolation factor on table
    f = (E - urr % energy(i_energy)) / &
         (urr % energy(i_energy + 1) - urr % energy(i_energy))

    ! sample probability table using the cumulative distribution
    r = prn()
    i_table = 1
    do
      if (urr % prob(i_energy, URR_CUM_PROB, i_table) > r) exit
      i_table = i_table + 1
    end do

    ! determine elastic, fission, and capture cross sections from probability
    ! table
    if (urr % interp == LINEAR_LINEAR) then
      elastic = (ONE - f) * urr % prob(i_energy, URR_ELASTIC, i_table) + &
           f * urr % prob(i_energy + 1, URR_ELASTIC, i_table)
      fission = (ONE - f) * urr % prob(i_energy, URR_FISSION, i_table) + &
           f * urr % prob(i_energy + 1, URR_FISSION, i_table)
      capture = (ONE - f) * urr % prob(i_energy, URR_N_GAMMA, i_table) + &
           f * urr % prob(i_energy + 1, URR_N_GAMMA, i_table)
    elseif (urr % interp == LOG_LOG) then
      ! Get logarithmic interpolation factor
      f = log(E / urr % energy(i_energy)) / &
           log(urr % energy(i_energy + 1) / urr % energy(i_energy))

      ! Calculate elastic cross section/factor
      elastic = ZERO
      if (urr % prob(i_energy, URR_ELASTIC, i_table) > ZERO) then
        elastic = exp((ONE - f) * log(urr % prob(i_energy, URR_ELASTIC, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_ELASTIC, &
             i_table)))
      end if

      ! Calculate fission cross section/factor
      fission = ZERO
      if (urr % prob(i_energy, URR_FISSION, i_table) > ZERO) then
        fission = exp((ONE - f) * log(urr % prob(i_energy, URR_FISSION, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_FISSION, &
             i_table)))
      end if

      ! Calculate capture cross section/factor
      capture = ZERO
      if (urr % prob(i_energy, URR_N_GAMMA, i_table) > ZERO) then
        capture = exp((ONE - f) * log(urr % prob(i_energy, URR_N_GAMMA, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_N_GAMMA, &
             i_table)))
      end if
    end if

    ! Determine treatment of inelastic scattering
    inelastic = ZERO
    if (urr % inelastic_flag > 0) then
      ! Get pointer to inelastic scattering reaction
      rxn => nuc % reactions(nuc % urr_inelastic)

      ! Get index on energy grid and interpolation factor
      i_energy = micro_xs(i_nuclide) % index_grid
      f = micro_xs(i_nuclide) % interp_factor

      ! Determine inelastic scattering cross section
      if (i_energy >= rxn % threshold) then
        inelastic = (ONE - f) * rxn % sigma(i_energy - rxn%threshold + 1) + &
             f * rxn % sigma(i_energy - rxn%threshold + 2)
      end if
    end if

    ! Multiply by smooth cross-section if needed
    if (urr % multiply_smooth) then
      elastic = elastic * micro_xs(i_nuclide) % elastic
      capture = capture * (micro_xs(i_nuclide) % absorption - &
           micro_xs(i_nuclide) % fission)
      fission = fission * micro_xs(i_nuclide) % fission
    end if

    ! Set elastic, absorption, fission, and total cross sections. Note that the
    ! total cross section is calculated as sum of partials rather than using the
    ! table-provided value
    micro_xs(i_nuclide) % elastic = elastic
    micro_xs(i_nuclide) % absorption = capture + fission
    micro_xs(i_nuclide) % fission = fission
    micro_xs(i_nuclide) % total = elastic + inelastic + capture + fission

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      micro_xs(i_nuclide) % nu_fission = nu_total(nuc, E) * &
           micro_xs(i_nuclide) % fission
    end if

  end subroutine calculate_urr_xs

!===============================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid at a certain
! energy
!===============================================================================

  function find_energy_index(E) result(e_index)

    real(8), intent(in) :: E ! energy of particle
    integer :: e_index ! unionized energy index

    ! if particle's energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (E < e_grid(1)) then
      e_index = 1
    elseif (E > e_grid(n_grid)) then
      e_index = n_grid - 1
    else
      e_index = binary_search(e_grid, n_grid, E)
    end if

  end function find_energy_index

end module cross_section
