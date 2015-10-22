module fdtd_data_cuda_module

use cudafor

use fdtd_data_module


implicit none

type fdtd_params_cuda
    integer, device, allocatable                 :: nx, ny, nz
    integer, device, allocatable                 :: nsrc
    integer, device, dimension(:,:), allocatable :: src
    real, device, allocatable                    :: dt !Length of the time step
    real, device, allocatable                    :: dx, dy, dz !Distance between 2 cells
    real, device, allocatable                    :: mu_0 !mu0, permeability of free space (in henry/meter)
    real, device, allocatable                    :: eps_0 !Epsilon0, permittivity of free space (in farad/meter)
    real, device, dimension(:), allocatable      :: jz
end type


type fdtd_field_cuda
    !FDTD field
    real, device, dimension(:,:,:), allocatable :: ex1, ex2, ex3
    real, device, dimension(:,:,:), allocatable :: ey1, ey2, ey3
    real, device, dimension(:,:,:), allocatable :: ez1, ez2, ez3

    real, device, dimension(:,:,:), allocatable :: hx
    real, device, dimension(:,:,:), allocatable :: hy
    real, device, dimension(:,:,:), allocatable :: hz

    real, device, dimension(:,:,:), allocatable :: dx1, dx2, dx3
    real, device, dimension(:,:,:), allocatable :: dy1, dy2, dy3
    real, device, dimension(:,:,:), allocatable :: dz1, dz2, dz3

    real, device, dimension(:,:,:), allocatable :: eps_i, eps_s
    real, device, dimension(:,:,:), allocatable :: tau_d, sigma

    !Mur boundary
    real, device, dimension(:,:,:), allocatable :: rp_x_1, rp_x_end
    real, device, dimension(:,:,:), allocatable :: rp_y_1, rp_y_end
    real, device, dimension(:,:,:), allocatable :: rp_z_1, rp_z_end
end type


contains

subroutine init_fdtd_parameters_cuda(params_cuda, params)
    !Input
    type(fdtd_params_cuda), allocatable, intent(inout) :: params_cuda
    type(fdtd_params), intent(in)                      :: params
    
    allocate(params_cuda)
    allocate(params_cuda%nx, params_cuda%ny, params_cuda%nz)
    allocate(params_cuda%nsrc)
    allocate(params_cuda%src(size(params%src, 1), size(params%src, 2)))
    allocate(params_cuda%dt, params_cuda%dx, params_cuda%dy, params_cuda%dz)
    allocate(params_cuda%mu_0, params_cuda%eps_0)
    allocate(params_cuda%jz(size(params%jz)))
    
    params_cuda%nx = params%nx
    params_cuda%ny = params%ny
    params_cuda%nz = params%nz
    params_cuda%nsrc = params%nsrc
    params_cuda%src = params%src
    params_cuda%dt = params%dt
    params_cuda%dx = params%dx
    params_cuda%dy = params%dy
    params_cuda%dz = params%dz
    params_cuda%mu_0 = params%mu_0
    params_cuda%eps_0 = params%eps_0
    params_cuda%jz = params%jz
end subroutine


subroutine init_fdtd_field_cuda(field_cuda, field)
    !Input
    type(fdtd_field_cuda), allocatable, intent(inout) :: field_cuda
    type(fdtd_field), intent(in)                      :: field
    
    allocate(field_cuda)

    allocate(field_cuda%ex1(size(field%ex1, 1), size(field%ex1, 2), size(field%ex1, 3)))
    allocate(field_cuda%ex2(size(field%ex2, 1), size(field%ex2, 2), size(field%ex2, 3)))
    allocate(field_cuda%ex3(size(field%ex3, 1), size(field%ex3, 2), size(field%ex3, 3)))
    allocate(field_cuda%ey1(size(field%ey1, 1), size(field%ey1, 2), size(field%ey1, 3)))
    allocate(field_cuda%ey2(size(field%ey2, 1), size(field%ey2, 2), size(field%ey2, 3)))
    allocate(field_cuda%ey3(size(field%ey3, 1), size(field%ey3, 2), size(field%ey3, 3)))
    allocate(field_cuda%ez1(size(field%ez1, 1), size(field%ez1, 2), size(field%ez1, 3)))
    allocate(field_cuda%ez2(size(field%ez2, 1), size(field%ez2, 2), size(field%ez2, 3)))
    allocate(field_cuda%ez3(size(field%ez3, 1), size(field%ez3, 2), size(field%ez3, 3)))

    allocate(field_cuda%hx(size(field%hx, 1), size(field%hx, 2), size(field%hx, 3)))
    allocate(field_cuda%hy(size(field%hy, 1), size(field%hy, 2), size(field%hy, 3)))
    allocate(field_cuda%hz(size(field%hz, 1), size(field%hz, 2), size(field%hz, 3)))

    allocate(field_cuda%dx1(size(field%dx1, 1), size(field%dx1, 2), size(field%dx1, 3)))
    allocate(field_cuda%dx2(size(field%dx2, 1), size(field%dx2, 2), size(field%dx2, 3)))
    allocate(field_cuda%dx3(size(field%dx3, 1), size(field%dx3, 2), size(field%dx3, 3)))
    allocate(field_cuda%dy1(size(field%dy1, 1), size(field%dy1, 2), size(field%dy1, 3)))
    allocate(field_cuda%dy2(size(field%dy2, 1), size(field%dy2, 2), size(field%dy2, 3)))
    allocate(field_cuda%dy3(size(field%dy3, 1), size(field%dy3, 2), size(field%dy3, 3)))
    allocate(field_cuda%dz1(size(field%dz1, 1), size(field%dz1, 2), size(field%dz1, 3)))
    allocate(field_cuda%dz2(size(field%dz2, 1), size(field%dz2, 2), size(field%dz2, 3)))
    allocate(field_cuda%dz3(size(field%dz3, 1), size(field%dz3, 2), size(field%dz3, 3)))

    allocate(field_cuda%eps_i(size(field%eps_i, 1), size(field%eps_i, 2), size(field%eps_i, 3)))
    allocate(field_cuda%eps_s(size(field%eps_s, 1), size(field%eps_s, 2), size(field%eps_s, 3)))
    allocate(field_cuda%tau_d(size(field%tau_d, 1), size(field%tau_d, 2), size(field%tau_d, 3)))
    allocate(field_cuda%sigma(size(field%sigma, 1), size(field%sigma, 2), size(field%sigma, 3)))

    allocate(field_cuda%rp_x_1(size(field%rp_x_1, 1), size(field%rp_x_1, 2), size(field%rp_x_1, 3)))
    allocate(field_cuda%rp_x_end(size(field%rp_x_end, 1), size(field%rp_x_end, 2), size(field%rp_x_end, 3)))
    allocate(field_cuda%rp_y_1(size(field%rp_y_1, 1), size(field%rp_y_1, 2), size(field%rp_y_1, 3)))
    allocate(field_cuda%rp_y_end(size(field%rp_y_end, 1), size(field%rp_y_end, 2), size(field%rp_y_end, 3)))
    allocate(field_cuda%rp_z_1(size(field%rp_z_1, 1), size(field%rp_z_1, 2), size(field%rp_z_1, 3)))
    allocate(field_cuda%rp_z_end(size(field%rp_z_end, 1), size(field%rp_z_end, 2), size(field%rp_z_end, 3)))

    field_cuda%ex1 = field%ex1
    field_cuda%ex2 = field%ex2
    field_cuda%ex3 = field%ex3
    field_cuda%ey1 = field%ey1
    field_cuda%ey2 = field%ey2
    field_cuda%ey3 = field%ey3
    field_cuda%ez1 = field%ez1
    field_cuda%ez2 = field%ez2
    field_cuda%ez3 = field%ez3

    field_cuda%hx = field%hx
    field_cuda%hy = field%hy
    field_cuda%hz = field%hz

    field_cuda%dx1 = field%dx1
    field_cuda%dx2 = field%dx2
    field_cuda%dx3 = field%dx3
    field_cuda%dy1 = field%dy1
    field_cuda%dy2 = field%dy2
    field_cuda%dy3 = field%dy3
    field_cuda%dz1 = field%dz1
    field_cuda%dz2 = field%dz2
    field_cuda%dz3 = field%dz3

    field_cuda%eps_i = field%eps_i
    field_cuda%eps_s = field%eps_s
    field_cuda%tau_d = field%tau_d
    field_cuda%sigma = field%sigma

    field_cuda%rp_x_1 = field%rp_x_1
    field_cuda%rp_x_end = field%rp_x_end
    field_cuda%rp_y_1 = field%rp_y_1
    field_cuda%rp_y_end = field%rp_y_end
    field_cuda%rp_z_1 = field%rp_z_1
    field_cuda%rp_z_end = field%rp_z_end
end subroutine


subroutine write_result_cuda(params, field, field_cuda, run_num, runs_count, output_path)
    !Input
    type(fdtd_params), intent(in)     :: params
    type(fdtd_field), intent(inout)   :: field
    type(fdtd_field_cuda), intent(in) :: field_cuda
    integer, intent(in)               :: run_num
    integer, intent(in)               :: runs_count
    character(len=*), intent(in)      :: output_path
    
    !Copy over results from GPU to RAM
    field%ex1 = field_cuda%ex1
    field%ex2 = field_cuda%ex2
    field%ex3 = field_cuda%ex3
    field%ey1 = field_cuda%ey1
    field%ey2 = field_cuda%ey2
    field%ey3 = field_cuda%ey3
    field%ez1 = field_cuda%ez1
    field%ez2 = field_cuda%ez2
    field%ez3 = field_cuda%ez3

    field%hx = field_cuda%hx
    field%hy = field_cuda%hy
    field%hz = field_cuda%hz
             
    field%dx1 = field_cuda%dx1
    field%dx2 = field_cuda%dx2
    field%dx3 = field_cuda%dx3
    field%dy1 = field_cuda%dy1
    field%dy2 = field_cuda%dy2
    field%dy3 = field_cuda%dy3
    field%dz1 = field_cuda%dz1
    field%dz2 = field_cuda%dz2
    field%dz3 = field_cuda%dz3

    field%eps_i = field_cuda%eps_i
    field%eps_s = field_cuda%eps_s
    field%tau_d = field_cuda%tau_d
    field%sigma = field_cuda%sigma

    field%rp_x_1 = field_cuda%rp_x_1
    field%rp_x_end = field_cuda%rp_x_end
    field%rp_y_1 = field_cuda%rp_y_1
    field%rp_y_end = field_cuda%rp_y_end
    field%rp_z_1 = field_cuda%rp_z_1
    field%rp_z_end = field_cuda%rp_z_end

    call write_result(params, field, run_num, runs_count, output_path)
end subroutine

end module
