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
    integer :: i, runs_count
    integer :: args_count
    character(len=5) :: is_cuda_arg !used to check if shall be ran using cuda
    logical :: is_cuda

    !Calculations data
    type(fdtd_params), allocatable :: params
    type(fdtd_field), allocatable  :: field
    
    !CUDA calculations data
    type(dim3)                          :: grid_size, block_size
    type(fdtd_params_cuda), allocatable :: params_cuda
    type(fdtd_field_cuda), allocatable  :: field_cuda
    integer :: cuda_stat
    
    !CUDA streams data
    integer :: err
    integer :: h_stream, d_stream!, e_stream
    !type(cudaEvent) :: h_event, d_event, e_event

    is_cuda = .false.

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
        print *, "Initializing GPU parameters..."
        call init_fdtd_parameters_cuda(params_cuda, params)
        print *, ""

        print *, "Initializing GPU field..."
        call init_fdtd_field_cuda(field_cuda, field)
        print *, ""
        
        block_size = dim3(TILE_SIZE, TILE_SIZE, TILE_SIZE)
        grid_size = dim3((params%nx + TILE_SIZE - 1)/TILE_SIZE, &
                         (params%nx + TILE_SIZE - 1)/TILE_SIZE, &
                         (params%nx + TILE_SIZE - 1)/TILE_SIZE)
                         
        !Stream data
        err = cudaStreamCreate(h_stream)
        err = cudaStreamCreate(d_stream)
        !err = cudaStreamCreate(e_stream)

        !err = cudaEventCreate(h_event)
        !err = cudaEventCreate(d_event)
        !err = cudaEventCreate(e_event)
    end if

    if(is_cuda) then
        print *, "Computation using GPU"
        print "(A, 3I4)", "Block size (x, y, z):", block_size%x, block_size%y, block_size%z
        print "(A, 3I6)", "Grid size (x, y, z):", grid_size%x, grid_size%y, grid_size%z
        
        g_nx = params%nx
        g_ny = params%ny
        g_nz = params%nz
        
        g_dt = params%dt
        g_dx = params%dx
        g_dy = params%dy
        g_dz = params%dz
        
        g_mu_0 = params%mu_0
        g_eps_0 = params%eps_0
        g_nsrc = params%nsrc
    else
        print *, "Computation using CPU"
    end if
    print *, ""
    
    !Main loop
    do i=1, params%runs_count, 3
        runs_count = (i-1)/3 + 1
        !First run
        print *, "Running iteration " // trim(str(i)) // "..."
        print *, ""
        
        !CUDA mode
        if(is_cuda) then
            !H stream
            call update_h_field_cuda<<<grid_size, block_size, 0, h_stream>>>(field_cuda%hx, field_cuda%hy, field_cuda%hz,                    &
                                                                             field_cuda%ex3, field_cuda%ey3, field_cuda%ez3)

            !err = cudaEventRecord(h_event, h_stream)


            !D stream
            !err = cudaStreamWaitEvent(d_stream, h_event, 0)

            call update_d_field_cuda<<<grid_size, block_size, 0, d_stream>>>(field_cuda%dx1, field_cuda%dy1, field_cuda%dz1, &
                                                                             field_cuda%dx3, field_cuda%dy3, field_cuda%dz3, &
                                                                             field_cuda%hx, field_cuda%hy, field_cuda%hz)

            call update_source_cuda<<<grid_size, block_size, 0, d_stream>>>(field_cuda%dz1, field_cuda%dz3,                                 &
                                                                            field_cuda%hx, field_cuda%hy,                                   &
                                                                            params_cuda%src, params_cuda%jz,                                &
                                                                            runs_count)

            !err = cudaEventRecord(d_event, d_stream)

            
            !E stream
            !err = cudaStreamWaitEvent(e_stream, d_event, 0)

            !call update_e_field_cuda<<<grid_size, block_size, 0, e_stream>>>(field_cuda%ex1, field_cuda%ey1, field_cuda%ez1, &
            !                                                                 field_cuda%ex3, field_cuda%ey3, field_cuda%ez3, &
            !                                                                 field_cuda%ex2, field_cuda%ey2, field_cuda%ez2, &
            !                                                                 field_cuda%dx1, field_cuda%dy1, field_cuda%dz1, &
            !                                                                 field_cuda%dx3, field_cuda%dy3, field_cuda%dz3, &
            !                                                                 field_cuda%dx2, field_cuda%dy2, field_cuda%dz2, &
            !                                                                 field_cuda%eps_i, field_cuda%eps_s,             &
            !                                                                 field_cuda%tau_d, field_cuda%sigma)
            
            !call update_mur_boundary_cuda<<<grid_size, block_size, 0, e_stream>>>(field_cuda%ex1, field_cuda%ey1, field_cuda%ez1,                 &
            !                                                                      field_cuda%ex3, field_cuda%ey3, field_cuda%ez3,                 &
            !                                                                      field_cuda%rp_x_1, field_cuda%rp_x_end,                         &
            !                                                                      field_cuda%rp_y_1, field_cuda%rp_y_end,                         &
            !                                                                      field_cuda%rp_z_1, field_cuda%rp_z_end)
                                                                                  
            !err = cudaEventRecord(e_event, e_stream)

            err = cudaMemcpyAsync(field%hx, field_cuda%hx, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
            err = cudaMemcpyAsync(field%hy, field_cuda%hy, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
            err = cudaMemcpyAsync(field%hz, field_cuda%hz, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)

            err = cudaMemcpyAsync(field%dx1, field_cuda%dx1, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
            err = cudaMemcpyAsync(field%dy1, field_cuda%dy1, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
            err = cudaMemcpyAsync(field%dz1, field_cuda%dz1, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
            
            !err = cudaMemcpyAsync(field%ex1, field_cuda%ex1, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
            !err = cudaMemcpyAsync(field%ey1, field_cuda%ey1, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
            !err = cudaMemcpyAsync(field%ez1, field_cuda%ez1, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
           
            !Write results
            !err = cudaStreamSynchronize(h_stream)
            !err = cudaStreamSynchronize(d_stream)
            !err = cudaStreamSynchronize(e_stream)

            err = cudaDeviceSynchronize()
            
            !call write_result(params, field,                   &
            !                  field%ex1, field%ey1, field%ez1, &
            !                  field%dx1, field%dy1, field%dz1, &
            !                  1, i, trim(params%output_path))
        !CPU mode
        else
            call update_h_field(params, field,                   &
            					field%ex3, field%ey3, field%ez3, &
					            1)
					            
            call update_d_field(params, field,                   &
            					field%dx1, field%dy1, field%dz1, &
					            field%dx3, field%dy3, field%dz3, &
					            1)
            
            call update_source(params, field,        &
                               field%dz1, field%dz3, &
                               1, runs_count)
            
            call update_e_field(params, field,                   &
                                field%ex1, field%ey1, field%ez1, &
                                field%ex3, field%ey3, field%ez3, &
                                field%ex2, field%ey2, field%ez2, &
                                field%dx1, field%dy1, field%dz1, &
                                field%dx3, field%dy3, field%dz3, &
                                field%dx2, field%dy2, field%dz2, &
                                1)
            
            call update_mur_boundary(params, field,                   &
                                     field%ex1, field%ey1, field%ez1, &
                                     field%ex3, field%ey3, field%ez3, &
                                     1)
        
            call write_result(params, field,                   &
                              field%ex1, field%ey1, field%ez1, &
                              field%dx1, field%dy1, field%dz1, &
                              1, i, trim(params%output_path))
        end if
        
!        !Second run
!        print *, "Running iteration " // trim(str(i+1)) // "..."
!        print *, ""
!        
!        !CUDA mode
!        if(is_cuda) then
!            !H stream
!            call update_h_field_cuda<<<grid_size, block_size, 0, h_stream>>>(field_cuda%hx, field_cuda%hy, field_cuda%hz,                    &
!                                                                             field_cuda%ex1, field_cuda%ey1, field_cuda%ez1)
!
!            err = cudaEventRecord(h_event, h_stream)
!
!            
!            !D stream
!            err = cudaStreamWaitEvent(d_stream, h_event, 0)
!
!            call update_d_field_cuda<<<grid_size, block_size, 0, d_stream>>>(field_cuda%dx2, field_cuda%dy2, field_cuda%dz2, &
!                                                                             field_cuda%dx1, field_cuda%dy1, field_cuda%dz1, &
!                                                                             field_cuda%hx, field_cuda%hy, field_cuda%hz)
!
!            call update_source_cuda<<<grid_size, block_size, 0, d_stream>>>(field_cuda%dz2, field_cuda%dz1,                                 &
!                                                                            field_cuda%hx, field_cuda%hy,                                   &
!                                                                            params_cuda%src, params_cuda%jz,                                &
!                                                                            runs_count)
!
!            err = cudaEventRecord(d_event, d_stream)
!
!
!            !E stream
!            err = cudaStreamWaitEvent(e_stream, d_event, 0)
!
!            call update_e_field_cuda<<<grid_size, block_size, 0, e_stream>>>(field_cuda%ex2, field_cuda%ey2, field_cuda%ez2, &
!                                                                             field_cuda%ex1, field_cuda%ey1, field_cuda%ez1, &
!                                                                             field_cuda%ex3, field_cuda%ey3, field_cuda%ez3, &
!                                                                             field_cuda%dx2, field_cuda%dy2, field_cuda%dz2, &
!                                                                             field_cuda%dx1, field_cuda%dy1, field_cuda%dz1, &
!                                                                             field_cuda%dx3, field_cuda%dy3, field_cuda%dz3, &
!                                                                             field_cuda%eps_i, field_cuda%eps_s,             &
!                                                                             field_cuda%tau_d, field_cuda%sigma)
!
!            call update_mur_boundary_cuda<<<grid_size, block_size, 0, e_stream>>>(field_cuda%ex2, field_cuda%ey2, field_cuda%ez2,                 &
!                                                                                  field_cuda%ex1, field_cuda%ey1, field_cuda%ez1,                 &
!                                                                                  field_cuda%rp_x_1, field_cuda%rp_x_end,                         &
!                                                                                  field_cuda%rp_y_1, field_cuda%rp_y_end,                         &
!                                                                                  field_cuda%rp_z_1, field_cuda%rp_z_end)
!
!            err = cudaEventRecord(e_event, e_stream)
!
!            err = cudaMemcpyAsync(field%hx, field_cuda%hx, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
!            err = cudaMemcpyAsync(field%hy, field_cuda%hy, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
!            err = cudaMemcpyAsync(field%hz, field_cuda%hz, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
!
!            err = cudaMemcpyAsync(field%dx2, field_cuda%dx2, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
!            err = cudaMemcpyAsync(field%dy2, field_cuda%dy2, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
!            err = cudaMemcpyAsync(field%dz2, field_cuda%dz2, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
!            
!            err = cudaMemcpyAsync(field%ex2, field_cuda%ex2, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
!            err = cudaMemcpyAsync(field%ey2, field_cuda%ey2, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
!            err = cudaMemcpyAsync(field%ez2, field_cuda%ez2, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
!
!            err = cudaStreamSynchronize(h_stream)
!            err = cudaStreamSynchronize(d_stream)
!            err = cudaStreamSynchronize(e_stream)
!            !err = cudaDeviceSynchronize()
!            
!            !Write results
!			call write_result(params, field,                   &
!                              field%ex2, field%ey2, field%ez2, &
!                              field%dx2, field%dy2, field%dz2, &
!                              2, i+1, trim(params%output_path))
!        !CPU mode
!        else
!            call update_h_field(params, field,                   &
!            					field%ex1, field%ey1, field%ez1, &
!            					2)
!
!            call update_d_field(params, field,                   &
!                        		field%dx2, field%dy2, field%dz2, &
!					            field%dx1, field%dy1, field%dz1, &
!					            2)
!
!            call update_source(params, field,        &
!                               field%dz2, field%dz1, &
!                               2, runs_count)
!            
!            call update_e_field(params, field,                   &
!                                field%ex2, field%ey2, field%ez2, &
!                                field%ex1, field%ey1, field%ez1, &
!                                field%ex3, field%ey3, field%ez3, &
!                                field%dx2, field%dy2, field%dz2, &
!                                field%dx1, field%dy1, field%dz1, &
!                                field%dx3, field%dy3, field%dz3, &
!                                2)
!            
!            call update_mur_boundary(params, field,                   &
!                                     field%ex2, field%ey2, field%ez2, &
!                                     field%ex1, field%ey1, field%ez1, &
!                                     2)
!
!			call write_result(params, field,                   &
!                              field%ex2, field%ey2, field%ez2, &
!                              field%dx2, field%dy2, field%dz2, &
!                              2, i+1, trim(params%output_path))
!        end if
!        
!        !Third run
!        print *, "Running iteration " // trim(str(i+2)) // "..."
!        print *, ""
!        
!        !CUDA mode
!        if(is_cuda) then
!            !H stream
!            call update_h_field_cuda<<<grid_size, block_size, 0, h_stream>>>(field_cuda%hx, field_cuda%hy, field_cuda%hz,                    &
!                                                                             field_cuda%ex2, field_cuda%ey2, field_cuda%ez2)
!
!            err = cudaEventRecord(h_event, h_stream)
!
!
!            !D stream
!            err = cudaStreamWaitEvent(d_stream, h_event, 0)
!
!            call update_d_field_cuda<<<grid_size, block_size, 0, d_stream>>>(field_cuda%dx3, field_cuda%dy3, field_cuda%dz3, &
!                                                                             field_cuda%dx2, field_cuda%dy2, field_cuda%dz2, &
!                                                                             field_cuda%hx, field_cuda%hy, field_cuda%hz)
!            
!            call update_source_cuda<<<grid_size, block_size, 0, d_stream>>>(field_cuda%dz3, field_cuda%dz2,                                 &
!                                                                            field_cuda%hx, field_cuda%hy,                                   &
!                                                                            params_cuda%src, params_cuda%jz,                                &
!                                                                            runs_count)
!
!            err = cudaEventRecord(d_event, d_stream)
!
!
!            !E stream
!            err = cudaStreamWaitEvent(e_stream, d_event, 0)
!
!            call update_e_field_cuda<<<grid_size, block_size, 0, e_stream>>>(field_cuda%ex3, field_cuda%ey3, field_cuda%ez3, &
!                                                                             field_cuda%ex2, field_cuda%ey2, field_cuda%ez2, &
!                                                                             field_cuda%ex1, field_cuda%ey1, field_cuda%ez1, &
!                                                                             field_cuda%dx3, field_cuda%dy3, field_cuda%dz3, &
!                                                                             field_cuda%dx2, field_cuda%dy2, field_cuda%dz2, &
!                                                                             field_cuda%dx1, field_cuda%dy1, field_cuda%dz1, &
!                                                                             field_cuda%eps_i, field_cuda%eps_s,             &
!                                                                             field_cuda%tau_d, field_cuda%sigma)
!
!            call update_mur_boundary_cuda<<<grid_size, block_size, 0, e_stream>>>(field_cuda%ex3, field_cuda%ey3, field_cuda%ez3,                 &
!                                                                                  field_cuda%ex2, field_cuda%ey2, field_cuda%ez1,                 &
!                                                                                  field_cuda%rp_x_1, field_cuda%rp_x_end,                         &
!                                                                                  field_cuda%rp_y_1, field_cuda%rp_y_end,                         &
!                                                                                  field_cuda%rp_z_1, field_cuda%rp_z_end)
!
!            err = cudaEventRecord(e_event, e_stream)
!
!            err = cudaMemcpyAsync(field%hx, field_cuda%hx, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
!            err = cudaMemcpyAsync(field%hy, field_cuda%hy, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
!            err = cudaMemcpyAsync(field%hz, field_cuda%hz, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, h_stream)
!
!            err = cudaMemcpyAsync(field%dx3, field_cuda%dx3, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
!            err = cudaMemcpyAsync(field%dy3, field_cuda%dy3, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
!            err = cudaMemcpyAsync(field%dz3, field_cuda%dz3, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, d_stream)
!           
!            err = cudaMemcpyAsync(field%ex3, field_cuda%ex3, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
!            err = cudaMemcpyAsync(field%ey3, field_cuda%ey3, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
!            err = cudaMemcpyAsync(field%ez3, field_cuda%ez3, params%nx*params%ny*params%nz, cudaMemcpyDeviceToHost, e_stream)
!            
!            !Write results
!            err = cudaStreamSynchronize(h_stream)
!            err = cudaStreamSynchronize(d_stream)
!            err = cudaStreamSynchronize(e_stream)
!            
!            !err = cudaDeviceSynchronize()
!
!			call write_result(params, field,                   &
!                              field%ex3, field%ey3, field%ez3, &
!                              field%dx3, field%dy3, field%dz3, &
!                              3, i+2, trim(params%output_path))
!        !CPU mode
!        else
!            call update_h_field(params, field,                   &
!            				    field%ex2, field%ey2, field%ez2, &
!            					3)
!
!            call update_d_field(params, field,                   &
!                                field%dx3, field%dy3, field%dz3, &
!					            field%dx2, field%dy2, field%dz2, &
!					            3)
!
!            call update_source(params, field,        &
!                               field%dz3, field%dz2, &
!                               3, runs_count)
!
!            call update_e_field(params, field,                   &
!                                field%ex3, field%ey3, field%ez3, &
!                                field%ex2, field%ey2, field%ez2, &
!                                field%ex1, field%ey1, field%ez1, &
!                                field%dx3, field%dy3, field%dz3, &
!                                field%dx2, field%dy2, field%dz2, &
!                                field%dx1, field%dy1, field%dz1, &
!                                3)
!
!            call update_mur_boundary(params, field,                   &
!                                     field%ex3, field%ey3, field%ez3, &
!                                     field%ex2, field%ey2, field%ez2, &
!                                     3)
!        
!			call write_result(params, field,                   &
!                              field%ex3, field%ey3, field%ez3, &
!                              field%dx3, field%dy3, field%dz3, &
!                              3, i+2, trim(params%output_path))
!        end if
    end do
end program
