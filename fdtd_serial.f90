include "utils_module.f90"
include "fdtd_data_module.f90"
include "fdtd_calculations_module.f90"


program fdtd_serial
use fdtd_data_module
use fdtd_calculations_module

implicit none
    !local vars
    integer :: i
    integer :: runs_count

    !calculations data
    type(fdtd_state), pointer :: state
    type(fdtd_field), pointer :: field

    !initialise data
    call init_fdtd_state(state, "data/input_params")
    call init_fdtd_field(field, state)
    call load_materials(state, field, "data/mat_specs_riken", "data/pgmdata/")
    
    !main loop
    !runs_count = int(real(state%t_max+2)/3.0)
    !do i=1, runs_count
        !first run
    !    call update_h_field(state, field, 0)
    !    call update_d_field(state, field, 0)
    !    call update_source(state, field, 0)
    !    call update_e_field(state, field, 0)
    !    call update_mur_boundary(state, field, 0)
        
        !second run
    !    call update_h_field(state, field, 1)
    !    call update_d_field(state, field, 1)
    !    call update_source(state, field, 1)
    !    call update_e_field(state, field, 1)
    !    call update_mur_boundary(state, field, 1)
        
        !third run
    !    call update_h_field(state, field, 2)
    !    call update_d_field(state, field, 2)
    !    call update_source(state, field, 2)
    !    call update_e_field(state, field, 2)
    !    call update_mur_boundary(state, field, 2)
    !end do
end program