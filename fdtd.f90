include "utils_module.f90"

include "fdtd_data_module.f90"
include "fdtd_calculations_module.f90"

include "fdtd_data_cuda_module.f90"
include "fdtd_calculations_cuda_module.f90"


program fdtd

use cudafor

use fdtd_data_module
use fdtd_calculations_module

use fdtd_data_cuda_module
use fdtd_calculations_cuda_module


implicit none

    !Local vars
    integer :: i
    integer :: args_count
    character(len=5) :: is_cuda_arg !used to check if shall be ran using cuda
    logical :: is_cuda

    !Calculations data
    type(fdtd_params), pointer :: params
    type(fdtd_field),  pointer :: field
    
    !CUDA calculations data
    type(dim3) :: grid_size, block_size
    device, type(fdtd_params_cuda) :: params_dev
    device, type(fdtd_field_cuda)  :: field_dev
    
    !Check if any command line arguments have been passed
    args_count = command_argument_count()
    if(args_count .ge. 1) then
        call get_command_argument(1, is_cuda_arg)
        
        !Make sure that there is only one
        if(args_count > 1 .or. is_cuda_arg /= "-cuda") then
            print *, "Usage: fdtd [-cuda]"
            stop
        end if
        
        is_cuda = .true.
    end if

    !Initialize data
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
    
    !Initialize CUDA
    if(is_cuda) then
        init_fdtd_parameters_cuda(params_dev, params)
        init_fdtd_field_cuda(field_dev, field)
        
        block_size = dim3(params%nx, 1, 1)
        grid_size = dim3(params%ny, params%nz, 1)
    end if
    
    !Main loop
    do i=1, params%runs_count, 3
        !First run
        print *, "Running iteration " // trim(str(i)) // "..."
        print *, ""
        
        if(is_cuda) then
            call update_h_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 1)
            call update_d_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 1)
            call update_source_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 1, (i-1)/3 + 1)
            call update_e_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 1)
            call update_mur_boundary_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 1)
            
            call write_result_cuda(params_dev, field_dev, 1, i, trim(params%output_path))
        else
            call update_h_field(params, field, 1)
            call update_d_field(params, field, 1)
            call update_source(params, field, 1, (i-1)/3 + 1)
            call update_e_field(params, field, 1)
            call update_mur_boundary(params, field, 1)
        
            call write_result(params, field, 1, i, trim(params%output_path))
        end if
        
        !Second run
        print *, "Running iteration " // trim(str(i+1)) // "..."
        print *, ""
        
        if(is_cuda) then
            call update_h_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 2)
            call update_d_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 2)
            call update_source_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 2, (i-1)/3 + 1)
            call update_e_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 2)
            call update_mur_boundary_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 2)
            
            call write_result_cuda(params_dev, field_dev, 2, i+1, trim(params%output_path))
        else
            call update_h_field(params, field, 2)
            call update_d_field(params, field, 2)
            call update_source(params, field, 2, (i-1)/3 + 1)
            call update_e_field(params, field, 2)
            call update_mur_boundary(params, field, 2)
        
            call write_result(params, field, 2, i+1, trim(params%output_path))
        end if
        
        !Third run
        print *, "Running iteration " // trim(str(i+2)) // "..."
        print *, ""
        
        if(is_cuda) then
            call update_h_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 3)
            call update_d_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 3)
            call update_source_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 3, (i-1)/3 + 1)
            call update_e_field_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 3)
            call update_mur_boundary_cuda<<<grid_size, block_size>>>(params_dev, field_dev, 3)
            
            call write_result_cuda(params_dev, field_dev, 3, i+2, trim(params%output_path))
        else
            call update_h_field(params, field, 3)
            call update_d_field(params, field, 3)
            call update_source(params, field, 3, (i-1)/3 + 1)
            call update_e_field(params, field, 3)
            call update_mur_boundary(params, field, 3)
        
            call write_result(params, field, 3, i+2, trim(params%output_path))
        end if
    end do

end program