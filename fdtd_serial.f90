include "utils_module.f90"
include "fdtd_data_module.f90"
include "fdtd_calculations_module.f90"


program fdtd_serial
use fdtd_data_module
use fdtd_calculations_module

implicit none
    !local vars
    integer :: i

    !calculations data
    type(fdtd_state), pointer :: state
    type(fdtd_field), pointer :: field

    !initialise data
    call init_fdtd_state(state, "data/input_params")
    call init_fdtd_field(field, state)
    call load_materials(state, field, "data/mat_specs_riken", trim(state%input_path))
    call setup_source(state, field)
    
    !main loop
    do i=1, state%runs_count, 3
        !first run
        call update_h_field(state, field, 1)
        call update_d_field(state, field, 1)
        call update_source(state, field, 1, (i-1)/3 + 1)
        call update_e_field(state, field, 1)
        call update_mur_boundary(state, field, 1)
        
        call write_result(state, field, 1, i, trim(state%output_path))
        
        !second run
        call update_h_field(state, field, 2)
        call update_d_field(state, field, 2)
        call update_source(state, field, 2, (i-1)/3 + 1)
        call update_e_field(state, field, 2)
        call update_mur_boundary(state, field, 2)
        
        call write_result(state, field, 2, i+1, trim(state%output_path))
        
        !third run
        call update_h_field(state, field, 3)
        call update_d_field(state, field, 3)
        call update_source(state, field, 3, (i-1)/3 + 1)
        call update_e_field(state, field, 3)
        call update_mur_boundary(state, field, 3)
        
        call write_result(state, field, 3, i+2, trim(state%output_path))
    end do
end program