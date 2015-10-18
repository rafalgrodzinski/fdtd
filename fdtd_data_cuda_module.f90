module fdtd_data_module
use fdtd_data_module


implicit none

type fdtd_params_cuda
    !Read from file
    integer                         :: nx, ny, nz
    integer                         :: runs_count
    character(len=128)              :: input_path
    character(len=128)              :: output_path
    integer                         :: elements_per_wave
    real                            :: wave_freq
    real                            :: pulse_width
    real                            :: pulse_modulation_freq
    integer                         :: nsrc
    device, integer, dimension(;,:) :: src
    real                            :: sigma
    real                            :: eps_s
    real                            :: eps_i
    real                            :: tau_d
    
    !Generated
    real :: pi = acos(-1.0) !Delicious pie
    real :: c = 3.0 * 10.0**8 !Light speed (v in m/s)
    real :: timeskip !Time step skip
    real :: lambda !Wave length (meters)
    real :: dt !Length of the time step
    real :: dx, dy, dz !Distance between 2 cells
    real :: mu_0 !mu0, permeability of free space (in henry/meter)
    real :: eps_0 !Epsilon0, permittivity of free space (in farad/meter)
    
    device, real, dimension(:) :: jz
end type


type fdtd_field_cuda
    !FDTD field
    device, real, dimension(:,:,:) :: ex1,ex2,ex3
    device, real, dimension(:,:,:) :: ey1,ey2,ey3 
    device, real, dimension(:,:,:) :: ez1,ez2,ez3 

    device, real, dimension(:,:,:) :: hx
    device, real, dimension(:,:,:) :: hy
    device, real, dimension(:,:,:) :: hz

    device, real, dimension(:,:,:) :: dx1,dx2,dx3 
    device, real, dimension(:,:,:) :: dy1,dy2,dy3 
    device, real, dimension(:,:,:) :: dz1,dz2,dz3 

    device, real, dimension(:,:,:) :: eps_i, eps_s
    device, real, dimension(:,:,:) :: tau_d, sigma

    !Mur boundary
    device, real, dimension(:,:,:) :: rp_x_1, rp_x_end
    device, real, dimension(:,:,:) :: rp_y_1, rp_y_end
    device, real, dimension(:,:,:) :: rp_z_1, rp_z_end
end type


subroutine init_fdtd_parameters_cuda(params, params_host)
    !Input
    device, type(fdtd_params_cuda), intent(inout) :: params
    type(fdtd_params), pointer, intent(in)        :: params_host
    
    !Local vars
    type(fdtd_params_cuda), pointer :: temp_params
    
    allocate(temp_params)
    allocate(temp_params%src(size(params%src)))
    allocate(temp_params%jz(size(params%jz)))
    
    temp_params%nx = params%nx
    temp_params%ny = params%ny
    temp_params%nz = params%nz
    temp_params%runs_count = params%runs_count
    temp_params%input_path = params%input_path
    temp_params%output_path = params%output_path
    temp_params%elements_per_wave = params%elements_per_wave
    temp_params%wave_freq = params%wave_freq
    temp_params%pulse_width = params%pulse_width
    temp_params%pulse_modulation_freq = params%pulse_modulation_freq
    temp_params%nsrc = params%nsrc
    temp_params%src = params%src
    temp_params%sigma = params%sigma
    temp_params%eps_s = params%eps_s
    temp_params%eps_i = params%eps_i
    temp_params%tau_d = params%tau_d
    temp_params%pi = params%pi
    temp_params%c = params%c
    temp_params%timeskip = params%timeskip
    temp_params%lambda = params%lambda
    temp_params%dt = params%dt
    temp_params%dx = params%dx
    temp_params%dy = params%dy
    temp_params%dz = params%dz
    temp_params%mu_0 = params%mu_0
    temp_params%eps_0 = params%eps_0
    temp_params%jz = params%jz
    
    allocate(params)
    params = temp_params
    
    deallocate(params)
end subroutine


subroutine init_fdtd_field_cuda(field_dev, field)
    !Input
    device, type(fdtd_field_cuda), intent(inout) :: field_dev
    type(fdtd_field), intent(in), pointer        :: field
    
    !Local vars
    type(fdtd_field_cuda), pointer :: temp_field
    
    allocate(temp_field)

    allocate(temp_field%ex1(sizeof(field%ex1)))
    allocate(temp_field%ex2(sizeof(field%ex2)))
    allocate(temp_field%ex3(sizeof(field%ex3)))
    allocate(temp_field%ey1(sizeof(field%ey1)))
    allocate(temp_field%ey2(sizeof(field%ey2)))
    allocate(temp_field%ey3(sizeof(field%ey3)))
    allocate(temp_field%ez1(sizeof(field%ez1)))
    allocate(temp_field%ez2(sizeof(field%ez2)))
    allocate(temp_field%ez3(sizeof(field%ez3)))

    allocate(temp_field%hx(sizeof(field%hx)))
    allocate(temp_field%hy(sizeof(field%hy)))
    allocate(temp_field%hz(sizeof(field%hz)))

    allocate(temp_field%dx1(sizeof(field%dx1)))
    allocate(temp_field%dx2(sizeof(field%dx2)))
    allocate(temp_field%dx3(sizeof(field%dx3)))
    allocate(temp_field%dy1(sizeof(field%dy1)))
    allocate(temp_field%dy2(sizeof(field%dy2)))
    allocate(temp_field%dy3(sizeof(field%dy3)))
    allocate(temp_field%dz1(sizeof(field%dz1)))
    allocate(temp_field%dz2(sizeof(field%dz2)))
    allocate(temp_field%dz3(sizeof(field%dz3)))

    allocate(temp_field%eps_i(sizeof(field%eps_i)))
    allocate(temp_field%eps_s(sizeof(field%eps_s)))
    allocate(temp_field%tau_d(sizeof(field%tau_d)))
    allocate(temp_field%sigma(sizeof(field%sigma)))

    allocate(temp_field%rp_x_1(sizeof(field%rp_x_1)))
    allocate(temp_field%rp_x_end(sizeof(field%rp_x_end)))
    allocate(temp_field%rp_y_1(sizeof(field%rp_y_1)))
    allocate(temp_field%rp_y_end(sizeof(field%rp_y_end)))
    allocate(temp_field%rp_z_1(sizeof(field%rp_z_1)))
    allocate(temp_field%rp_z_end(sizeof(field%rp_z_end)))

    temp_field%ex1 = field%ex1
    temp_field%ex2 = field%ex2
    temp_field%ex3 = field%ex3
    temp_field%ey1 = field%ey1
    temp_field%ey2 = field%ey2
    temp_field%ey3 = field%ey3
    temp_field%ez1 = field%ez1
    temp_field%ez2 = field%ez2
    temp_field%ez3 = field%ez3
                  
    temp_field%hx = field%hx
    temp_field%hy = field%hy
    temp_field%hz = field%hz
                  
    temp_field%dx1 = field%dx1
    temp_field%dx2 = field%dx2
    temp_field%dx3 = field%dx3
    temp_field%dy1 = field%dy1
    temp_field%dy2 = field%dy2
    temp_field%dy3 = field%dy3
    temp_field%dz1 = field%dz1
    temp_field%dz2 = field%dz2
    temp_field%dz3 = field%dz3

    temp_field%eps_i = field%eps_i
    temp_field%eps_s = field%eps_s
    temp_field%tau_d = field%tau_d
    temp_field%sigma = field%sigma

    temp_field%rp_x_1 = field%rp_x_1
    temp_field%rp_x_end = field%rp_x_end
    temp_field%rp_y_1 = field%rp_y_1
    temp_field%rp_y_end = field%rp_y_end
    temp_field%rp_z_1 = field%rp_z_1
    temp_field%rp_z_end = field%rp_z_end
    
    allocate(field_dev)
    field_dev = temp_field
    
    deallocate(temp_field)
end subroutine


subroutine write_result_cuda(params, field, field_dev, run_num, runs_count, output_path)
    !Input
    type(fdtd_params), pointer, intent(in)    :: params
    type(fdtd_field), pointer, intent(in)     :: field
    device, type(fdtd_field_cuda), intent(in) :: field_dev
    integer, intent(in)                       :: run_num
    integer, intent(in)                       :: runs_count
    character(len=*), intent(in)              :: output_path
    
    !Local vars
    type(fdtd_field_cuda), pointer :: temp_field    

    !Copy over results from GPU to CPU
    allocate(temp_field)
    temp_field = field_dev
    
    field%ex1 = temp_field%ex1
    field%ex2 = temp_field%ex2
    field%ex3 = temp_field%ex3
    field%ey1 = temp_field%ey1
    field%ey2 = temp_field%ey2
    field%ey3 = temp_field%ey3
    field%ez1 = temp_field%ez1
    field%ez2 = temp_field%ez2
    field%ez3 = temp_field%ez3
                  
    field%hx = temp_field%hx
    field%hy = temp_field%hy
    field%hz = temp_field%hz
             
    field%dx1 = temp_field%dx1
    field%dx2 = temp_field%dx2
    field%dx3 = temp_field%dx3
    field%dy1 = temp_field%dy1
    field%dy2 = temp_field%dy2
    field%dy3 = temp_field%dy3
    field%dz1 = temp_field%dz1
    field%dz2 = temp_field%dz2
    field%dz3 = temp_field%dz3

    field%eps_i = temp_field%eps_i
    field%eps_s = temp_field%eps_s
    field%tau_d = temp_field%tau_d
    field%sigma = temp_field%sigma

    field%rp_x_1 = temp_field%rp_x_1
    field%rp_x_end = temp_field%rp_x_end
    field%rp_y_1 = temp_field%rp_y_1
    field%rp_y_end = temp_field%rp_y_end
    field%rp_z_1 = temp_field%rp_z_1
    field%rp_z_end = temp_field%rp_z_end

    deallocate(temp_field)
    
    write_result(params, field, run_num, runs_count, output_path)
end subroutine

end module
