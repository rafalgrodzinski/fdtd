include "utils_module.f90"

include "fdtd_data_module.f90"
include "fdtd_calculations_module.f90"

include "fdtd_data_cuda_module.f90"
include "fdtd_calculations_cuda_module.f90"


program fdtd
use fdtd_data_module
use fdtd_calculations_module

use fdtd_data_cuda_module
use fdtd_calculations_cuda_module

implicit none
    !local vars
    integer :: i
    integer :: args_count
    character(len=5) :: is_cuda_arg !used to check if shall be ran using cuda
    logical :: is_cuda
    

    !calculations data
    type(fdtd_params), pointer :: params
    type(fdtd_field), pointer  :: field
    
    !check if any command line arguments have been passed
    args_count = command_argument_count()
    if(args_count .ge. 1) then
        call get_command_argument(1, is_cuda_arg)
        
        !make sure that there is only one
        if(args_count > 1 .or. is_cuda_arg /= "-cuda") then
            print *, "Usage: fdtd [-cuda]"
            stop
        end if
        
        is_cuda = .true.
    end if

    !initialise data
    print *, "Reading parameters..."
    call init_fdtd_parameters(params, "data/input_params")
    call print_parameters(params)
    print *, ""
    
    print *, "Initialising field..."
    call init_fdtd_field(field, params)
    print *, ""
    
    print *, "Reading pgm data..."
    call load_materials(params, field, "data/mat_specs_riken", trim(params%input_path))
    call setup_source(params, field)
    print *, ""
    
    !main loop
    do i=1, params%runs_count, 3
        !first run
        print *, "Running iteration " // trim(str(i)) // "..."
        print *, ""
        call update_h_field(params, field, 1)
        call update_d_field(params, field, 1)
        call update_source(params, field, 1, (i-1)/3 + 1)
        call update_e_field(params, field, 1)
        call update_mur_boundary(params, field, 1)
        
        call write_result(params, field, 1, i, trim(params%output_path))
        
        !second run
        print *, "Running iteration " // trim(str(i+1)) // "..."
        print *, ""
        call update_h_field(params, field, 2)
        call update_d_field(params, field, 2)
        call update_source(params, field, 2, (i-1)/3 + 1)
        call update_e_field(params, field, 2)
        call update_mur_boundary(params, field, 2)
        
        call write_result(params, field, 2, i+1, trim(params%output_path))
        
        !third run
        print *, "Running iteration " // trim(str(i+2)) // "..."
        print *, ""
        call update_h_field(params, field, 3)
        call update_d_field(params, field, 3)
        call update_source(params, field, 3, (i-1)/3 + 1)
        call update_e_field(params, field, 3)
        call update_mur_boundary(params, field, 3)
        
        call write_result(params, field, 3, i+2, trim(params%output_path))
    end do

end program