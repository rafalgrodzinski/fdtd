module fdtd_data_module

use utils_module


implicit none

type fdtd_params
    !Read from file
    integer              :: nx, ny, nz
    integer              :: runs_count
    character(len=128)   :: input_path
    character(len=128)   :: output_path
    integer              :: elements_per_wave
    real                 :: wave_freq
    real                 :: pulse_width
    real                 :: pulse_modulation_freq
    integer              :: nsrc
    integer, allocatable :: src(:,:)
    real                 :: sigma
    real                 :: eps_s
    real                 :: eps_i
    real                 :: tau_d
    
    !Generated
    real :: pi = acos(-1.0) !Delicious pie
    real :: c = 3.0 * 10.0**8 !Light speed (v in m/s)
    real :: timeskip !Time step skip
    real :: lambda !Wave length (meters)
    real :: dt !Length of the time step
    real :: dx, dy, dz !Distance between 2 cells
    real :: mu_0 !mu0, permeability of free space (in henry/meter)
    real :: eps_0 !Epsilon0, permittivity of free space (in farad/meter)
    
    real, dimension(:), allocatable :: jz
end type


type fdtd_field
    !FDTD field
    real, pinned, dimension(:,:,:), allocatable :: ex1, ex2, ex3
    real, pinned, dimension(:,:,:), allocatable :: ey1, ey2, ey3 
    real, pinned, dimension(:,:,:), allocatable :: ez1, ez2, ez3 

    real, pinned, dimension(:,:,:), allocatable :: hx, hy, hz

    real, pinned, dimension(:,:,:), allocatable :: dx1, dx2, dx3 
    real, pinned, dimension(:,:,:), allocatable :: dy1, dy2, dy3 
    real, pinned, dimension(:,:,:), allocatable :: dz1, dz2, dz3 

    real, dimension(:,:,:), allocatable :: eps_i, eps_s
    real, dimension(:,:,:), allocatable :: tau_d, sigma

    !Mur boundary
    real, dimension(:,:,:), allocatable :: rp_x_1, rp_x_end
    real, dimension(:,:,:), allocatable :: rp_y_1, rp_y_end
    real, dimension(:,:,:), allocatable :: rp_z_1, rp_z_end
end type


contains

subroutine init_fdtd_parameters(params, file_name)
    !Input
    type(fdtd_params), allocatable, intent(inout) :: params
    character(len=*), intent(in)                  :: file_name
    
    !Local vars
    integer, parameter :: file_unit = 51
    integer            :: error_code
    character          :: temp_c
    integer            :: temp_i
    integer            :: i
    
    allocate(params)
    
    open(unit=file_unit, file=file_name, status="old", iostat=error_code)
    call check_error(error_code, "Couldn't open file "//file_name)

    !nx_ny_nz (field size)
    read(file_unit, *) temp_c, params%nx, params%ny, params%nz
    !t_max (simulation runs count)
    read(file_unit, *) temp_c, temp_i
    params%runs_count = ((temp_i-1)/3 + 1)*3 !runs count has to be divisible by 3
    !unused (nf)
    read(file_unit, *) temp_c, temp_i
    !env_set_dir (input path)
    read(file_unit, *) temp_c, params%input_path
    !unused (env_file_prefix)
    read(file_unit, *) temp_c, temp_c
    !output_dir (output path)
    read(file_unit, *) temp_c, params%output_path
    !unused (output_format)
    read(file_unit, *) temp_c, temp_c
    !unused (impulse_resp_flag)
    read(file_unit, *) temp_c, temp_i
    !unused (pec_flag)
    read(file_unit, *) temp_c, temp_i
    !unused (read_env_flag)
    read(file_unit, *) temp_c, temp_i
    !unused (output_flag)
    read(file_unit, *) temp_c, temp_i
    !unused (bzip2_flag)
    read(file_unit, *) temp_c, temp_i
    !unused (output_start)
    read(file_unit, *) temp_c, temp_i, temp_i, temp_i
    !unused (output_finish)
    read(file_unit, *) temp_c, temp_i, temp_i, temp_i
    !unused (source_type)
    read(file_unit, *) temp_c, temp_c
    !elements_per_wavelength
    read(file_unit, *) temp_c, params%elements_per_wave
    !wave_freq
    read(file_unit, *) temp_c, params%wave_freq
    !pulse_width
    read(file_unit, *) temp_c, params%pulse_width
    !pulse_modulation_frequency
    read(file_unit, *) temp_c, params%pulse_modulation_freq
    !number_of_excitation_sources
    read(file_unit, *) temp_c, params%nsrc
    allocate (params%src(params%nsrc, 1:3))
    !source_location
    do i = 1, params%nsrc
      read(file_unit, *) temp_c, params%src(i, 1:3)
    enddo
    !unused (pulse_type)
    read(file_unit, *) temp_c, temp_i
    !fsigma (sigma)
    read(file_unit, *) temp_c, params%sigma
    !feps_s (eps_s)
    read(file_unit, *) temp_c, params%eps_s
    !feps_inf (eps_i)
    read(file_unit, *) temp_c, params%eps_i
    !ftau_d (tau_d)
    read(file_unit, *) temp_c, params%tau_d
    
    close(file_unit)
    
    !generate rest of the values
    params%timeskip = 1.0
    params%lambda = params%c / params%wave_freq
    params%dx = params%lambda/params%elements_per_wave
    params%dy = params%dx
    params%dz = params%dx
    params%dt = 1.0d0 * params%timeskip / (params%c * sqrt(1.0d0/(params%dx**2) + 1.0d0/(params%dy**2) + 1.0d0/(params%dz**2)))
    params%mu_0 = 4 * params%pi*10**(-7.0)
    params%eps_0 = 1.0/(params%mu_0 * params%c * params%c)
end subroutine


subroutine delete_fdtd_parameters(params)
    !Input
    type(fdtd_params), allocatable, intent(inout) :: params
    
    deallocate(params%src)
    deallocate(params%jz)
    deallocate(params)
end subroutine


subroutine init_fdtd_field(field, params)
    !Input
    type(fdtd_field), allocatable, intent(inout)  :: field
    type(fdtd_params), allocatable, intent(inout) :: params
    
    !Local vars
    integer :: ix, iy, iz
    
    allocate(field)
    
    !Initialise H field
    allocate(field%hx(params%nx, params%ny, params%nz))
    allocate(field%hy(params%nx, params%ny, params%nz))
    allocate(field%hz(params%nx, params%ny, params%nz))

    field%hx = 0.0
    field%hy = 0.0
    field%hz = 0.0
    
    !Initialise D field
    allocate(field%dx1(params%nx, params%ny, params%nz))
    allocate(field%dx2(params%nx, params%ny, params%nz))
    allocate(field%dx3(params%nx, params%ny, params%nz))
    
    field%dx1 = 0.0
    field%dx2 = 0.0
    field%dx3 = 0.0

    allocate(field%dy1(params%nx, params%ny, params%nz))
    allocate(field%dy2(params%nx, params%ny, params%nz))
    allocate(field%dy3(params%nx, params%ny, params%nz))
    
    field%dy1 = 0.0
    field%dy2 = 0.0
    field%dy3 = 0.0

    allocate(field%dz1(params%nx, params%ny, params%nz))
    allocate(field%dz2(params%nx, params%ny, params%nz))
    allocate(field%dz3(params%nx, params%ny, params%nz))
    
    field%dz1 = 0.0
    field%dz2 = 0.0
    field%dz3 = 0.0

    !Initialise E field
    allocate(field%ex1(params%nx, params%ny, params%nz))
    allocate(field%ex2(params%nx, params%ny, params%nz))
    allocate(field%ex3(params%nx, params%ny, params%nz))
    
    field%ex1 = 0.0
    field%ex2 = 0.0
    field%ex3 = 0.0

    allocate(field%ey1(params%nx, params%ny, params%nz))
    allocate(field%ey2(params%nx, params%ny, params%nz))
    allocate(field%ey3(params%nx, params%ny, params%nz))
    
    field%ey1 = 0.0
    field%ey2 = 0.0
    field%ey3 = 0.0

    allocate(field%ez1(params%nx, params%ny, params%nz))
    allocate(field%ez2(params%nx, params%ny, params%nz))
    allocate(field%ez3(params%nx, params%ny, params%nz))
    
    field%ez1 = 0.0
    field%ez2 = 0.0
    field%ez3 = 0.0

    allocate(field%eps_i(params%nx, params%ny, params%nz))
    allocate(field%eps_s(params%nx, params%ny, params%nz))

    allocate(field%tau_d(params%nx, params%ny, params%nz))
    allocate(field%sigma(params%nx, params%ny, params%nz))

    field%eps_i = 0.0
    field%eps_s = 0.0
    field%tau_d = 0.0
    field%sigma = 0.0
    
    field%sigma = params%sigma
    field%eps_s = params%eps_s
    field%eps_i = params%eps_i
    field%tau_d = params%tau_d
    
    !Initialise mur boundary
    allocate(field%rp_x_1(2, params%ny, params%nz))
    allocate(field%rp_x_end(params%nx-1:params%nx, params%ny, params%nz))
    
    allocate(field%rp_y_1(params%nx, 2, params%nz))
    allocate(field%rp_y_end(params%nx, params%ny-1:params%ny, params%nz))
    
    allocate(field%rp_z_1(params%nx, params%ny, 2))
    allocate(field%rp_z_end(params%nx, params%ny, params%nz-1:params%nz))
    
    field%rp_z_1 = 0.0
    field%rp_z_end = 0.0

    call setup_mur_boundary(params, field)
end subroutine


subroutine delete_fdtd_field(field)
    !Input
    type(fdtd_field), allocatable, intent(inout) :: field
    
    deallocate(field%ex1)
    deallocate(field%ex2)
    deallocate(field%ex3)

    deallocate(field%ey1)
    deallocate(field%ey2)
    deallocate(field%ey3)

    deallocate(field%ez1)
    deallocate(field%ez2)
    deallocate(field%ez3)

    deallocate(field%hx)
    deallocate(field%hy)
    deallocate(field%hz)

    deallocate(field%dx1)
    deallocate(field%dx2)
    deallocate(field%dx3)

    deallocate(field%dy1)
    deallocate(field%dy2)
    deallocate(field%dy3)

    deallocate(field%dz1)
    deallocate(field%dz2)
    deallocate(field%dz3)

    deallocate(field%eps_i)
    deallocate(field%eps_s)

    deallocate(field%tau_d)
    deallocate(field%sigma)
    
    deallocate(field%rp_x_1)
    deallocate(field%rp_x_end)

    deallocate(field%rp_y_1)
    deallocate(field%rp_y_end)

    deallocate(field%rp_z_1)
    deallocate(field%rp_z_end)
    
    deallocate(field)
end subroutine


subroutine load_materials(params, field, specs_file_name, materials_path)
    !Input
    type(fdtd_params), intent(inout) :: params
    type(fdtd_field), intent(inout)  :: field  
    character(len=*), intent(in)     :: specs_file_name
    character(len=*), intent(in)     :: materials_path
    
    !Local vars
    integer, parameter                :: file_unit = 100
    integer                           :: error_code
    real, dimension(:,:), allocatable :: specs
    integer, parameter                :: specs_count = 94
    integer                           :: spec_code
    character                         :: dummy_char
    integer                           :: ix, iy, iz
    character(len=128)                :: material_file_name
    integer                           :: material_width, material_height
    
    allocate(specs(0:specs_count, 1:4))
    specs = 0.0
    
    !load material specs
    open(file_unit, file=specs_file_name, status="old", iostat=error_code)
    call check_error(error_code, "Couldn't open file " // specs_file_name)
    
    spec_code = 0
    do while(spec_code .lt. specs_count)
        read(file_unit, *) spec_code, dummy_char, specs(spec_code, 1:4)
    end do
    
    close(unit=file_unit)
    
    !load materials
    do iz=1, params%nz
        !generte file name, starting with v1_00001.pgm
        material_file_name = generate_file_name(materials_path // "v1_", ".pgm", iz)
        open(unit=file_unit, file=material_file_name, status="old", iostat=error_code)
        call check_error(error_code, "Couldn't open file "//material_file_name)
        
        read(file_unit, *) dummy_char, dummy_char, dummy_char, material_width, material_height, dummy_char
        
        do iy=1, params%ny
            do ix=1, params%nx
                read(file_unit, *) spec_code
            
                field%sigma(ix, iy, iz) = specs(spec_code, 1)
                field%eps_s(ix, iy, iz) = specs(spec_code, 2)
                field%eps_i(ix, iy, iz) = specs(spec_code, 3)
                field%tau_d(ix, iy, iz) = specs(spec_code, 4)
            end do
        end do
        
        close(unit=file_unit)
    end do
    
    deallocate(specs)
end subroutine


subroutine setup_mur_boundary(params, field)
    !Input
    type(fdtd_params), intent(inout) :: params
    type(fdtd_field),  intent(inout) :: field
    
    !Local vars
    integer :: ix, iy, iz
    
    !Setup rp_x
    field%rp_x_1 = 0.0
    field%rp_x_end = 0.0

    do iz=1, params%nz
	    do iy=1, params%ny
            do ix=1, 2
                field%rp_x_1(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                               &
                                            (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                       &
                                            (1 + cmplx(0.0, 2 * params%pi * params%wave_freq * field%tau_d(ix, iy,iz))) - &
                                            cmplx(0, field%sigma(ix, iy, iz) /                                          &
                                             (2 * params%pi * params%wave_freq * params%eps_0)))
            end do
            
            do ix=params%nx-1, params%nx 
                field%rp_x_end(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                  &
                                                  (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                      &
                                                  (1 + cmplx(0, 2 * params%pi * params%wave_freq * field%tau_d(ix, iy, iz))) - &
                                                  cmplx(0, field%sigma(ix, iy, iz) /                                         &
                                                   (2 * params%pi * params%wave_freq * params%eps_0)))
            end do
	    end do
    end do
    
    !Setup rp_y
    field%rp_y_1 = 0.0
    field%rp_y_end = 0.0

    do iz=1, params%nz
        do ix=1, params%nx 
            do iy=1, 2 
                field%rp_y_1(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                      &
                                                (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                          &
                                                (1.0 + cmplx(0, 2.0 * params%pi * params%wave_freq * field%tau_d(ix, iy, iz))) - &
                                                cmplx(0, field%sigma(ix, iy, iz) /                                             &
                                                (2.0 * params%pi * params%wave_freq * params%eps_0)))
            end do
	    
            do iy=params%ny-1, params%ny
                field%rp_y_end(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                      &
                                                  (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                          &
                                                  (1.0 + cmplx(0, 2.0 * params%pi * params%wave_freq * field%tau_d(ix, iy, iz))) - &
                                                  cmplx(0, field%sigma(ix, iy, iz) /                                             &
                                                  (2 * params%pi * params%wave_freq * params%eps_0)))
            end do
	    end do
    end do
    
    !Setup rp_z
    field%rp_z_1 = 0.0
    field%rp_z_end = 0.0

    do iy=1, params%ny
        do ix=1, params%nx
            do iz=1, 2 
                field%rp_z_1(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                  &
                                                (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                      &
                                                (1 + cmplx(0, 2 * params%pi * params%wave_freq * field%tau_d(ix, iy, iz))) - &
                                                cmplx(0, field%sigma(ix, iy, iz) /                                         &
                                                (2 * params%pi * params%wave_freq * params%eps_0)))
            end do
            
            do iz=params%nz-1, params%nz 
                field%rp_z_end(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                  &
                                                  (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                      &
                                                  (1 + cmplx(0, 2 * params%pi * params%wave_freq * field%tau_d(ix, iy, iz))) - &
                                                  cmplx(0, field%sigma(ix, iy, iz) /                                         &
                                                  (2 * params%pi * params%wave_freq * params%eps_0)))
            end do
	    end do
    end do
end subroutine


subroutine setup_source(params, field)
    !Input
    type(fdtd_params), intent(inout) :: params
    type(fdtd_field), intent(inout)  :: field

    !Local vars
    integer :: fine
    integer :: temp
    integer :: i, i2
    integer :: istart
    real, dimension(:), allocatable :: tmpdata, tmpdata2
    
    allocate(params%jz(1:2**16))
    allocate(tmpdata(-2**16:2**16))
    allocate(tmpdata2(-2**16:2**16))
    
    !fine & temp
    fine = int(2**13 * params%pulse_width * params%wave_freq * params%dt)

    temp = 1/(params%pulse_width * params%wave_freq)/(params%dt / fine)/2
    
    !tmpdata
    do i=-2**14, 2**14
        tmpdata(i)=cmplx(0.0)
    end do
      
    do i=-temp-1,temp+1
         tmpdata(i)= exp(-(real(i)/((real(temp)+1.0)/4.0))**2)
    end do
      
    do i = -temp-1,temp+1
      tmpdata(i) = tmpdata(i) *cos(2.0*acos(-1.0)* params%pulse_modulation_freq * params%wave_freq*i*(params%dt/fine))     
    enddo

    !istart
    do i=-2**12,2**12-1
         if( (abs(tmpdata(i)).gt.1e-9).and.(mod(i,fine).eq.0)) then
            istart = i 
            exit
         endif
    enddo
    
    !setup jz 1/2
    params%jz = 0.0
    
    i2 = 0
    do i=istart, temp+1, fine
          i2=i2+1
          params%jz(i2)=tmpdata(i)*10**(-15.)/params%dt/3.0
    end do
    
    !setup tmpdata2
    do i=2, 2**14
        tmpdata2(i-1)=(real((params%jz(i+1)-params%jz(i))/params%dt)+real((params%jz(i)-params%jz(i-1))/params%dt)) &
        /2.0*(params%dt*params%dz)/(params%dx*params%dy*params%dz)
    end do
    
    !setup jz 2/2
    do i=1, 2**14
        params%jz(i) = tmpdata2(i)
    end do
end subroutine


subroutine print_parameters(params)
    !Input
    type(fdtd_params), intent(inout) :: params

    !Local vars
    integer :: i
    
    print "(A, 3I4)",   "Field size (x, y, z):      ", params%nx, params%ny, params%nz
    print "(A, I4)",    "Iterations count:          ", params%runs_count
    print "(A, A)",     "Input path:                ", trim(params%input_path)
    print "(A, A)",     "Output path:               ", trim(params%output_path)
    print "(A, I4)",    "Elements per wavelength:   ", params%elements_per_wave
    print "(A, E11.3)", "Wave frequency:            ", params%wave_freq
    print "(A, E11.3)", "Pulse width:               ", params%pulse_width
    print "(A, E11.3)", "Pulse modulation frequency:", params%pulse_modulation_freq
    do i=1, params%nsrc
        print "(A, 3I4)", "Source position (x, y, z): ", params%src(i, 1:3)
    end do
    print "(A, E11.3)", "Default sigma:             ", params%sigma
    print "(A, E11.3)", "Default eps_s:             ", params%eps_s
    print "(A, E11.3)", "Default eps_i:             ", params%eps_i
    print "(A, E11.3)", "Default tau_d:             ", params%tau_d
end subroutine


subroutine write_result(params, field,                   &
                        ex_source, ey_source, ez_source, &
                        dx_source, dy_source, dz_source, &
                        run_num, runs_count, output_path)

    !Input
    type(fdtd_params), intent(in)      :: params
    type(fdtd_field), intent(in)       :: field
    
	real, dimension(:,:,:), intent(in) :: ex_source
    real, dimension(:,:,:), intent(in) :: ey_source
    real, dimension(:,:,:), intent(in) :: ez_source
    
    real, dimension(:,:,:), intent(in) :: dx_source
    real, dimension(:,:,:), intent(in) :: dy_source
    real, dimension(:,:,:), intent(in) :: dz_source
    
    integer, intent(in)                :: run_num
    integer, intent(in)                :: runs_count
    character(len=*), intent(in)       :: output_path
    
    !Local vars
    integer, parameter :: file_unit = 100
    character(len=128) :: output_file_name
    integer            :: error_code
    integer            :: i
    integer            :: ix, iy, iz, count

    !Output x
    !Generte file name, starting with E_field_x_00001.out
    output_file_name = generate_file_name(output_path // "E_field_x_", ".out", runs_count)
    open(file_unit, file=output_file_name, status="new", access="sequential", form="formatted", iostat=error_code)
    call check_error(error_code, "Couldn't create file " // output_file_name)
    
    do i=1, params%nsrc
        iy = params%src(i, 2)
        iz = params%src(i, 3)
        do ix=1, params%nx
            write(file_unit, '(3I4, 9E11.3)') ix, iy, iz,                                                          &
                                              dx_source(ix, iy, iz), dy_source(ix, iy, iz), dz_source(ix, iy, iz), &
                                              field%hx(ix, iy, iz),  field%hy(ix, iy, iz),  field%hz(ix, iy, iz),  &
                                              ex_source(ix, iy, iz), ey_source(ix, iy, iz), ez_source(ix, iy, iz)
        end do
    end do
    
    close(file_unit)
    
    !Output y
    !generte file name, starting with E_field_y_00001.out
    output_file_name = generate_file_name(output_path // "E_field_y_", ".out", runs_count)
    open(file_unit, file=output_file_name, status="new", access="sequential", form="formatted", iostat=error_code)
    call check_error(error_code, "Couldn't create file " // output_file_name)
    
    do i=1, params%nsrc
        ix = params%src(i, 1)
        iz = params%src(i, 3)
        do iy=1, params%ny
            write(file_unit, '(3I4, 9E11.3)') ix, iy, iz,                                                          &
                                              dx_source(ix, iy, iz), dy_source(ix, iy, iz), dz_source(ix, iy, iz), &
                                              field%hx(ix, iy, iz),  field%hy(ix, iy, iz),  field%hz(ix, iy, iz),  &
                                              ex_source(ix, iy, iz), ey_source(ix, iy, iz), ez_source(ix, iy, iz)
        end do
    end do
    
    close(file_unit)
    
    !Output z
    !generte file name, starting with E_field_z_00001.out
    output_file_name = generate_file_name(output_path // "E_field_z_", ".out", runs_count)
    open(file_unit, file=output_file_name, status="new", access="sequential", form="formatted", iostat=error_code)
    call check_error(error_code, "Couldn't create file " // output_file_name)
    
    do i=1, params%nsrc
        ix = params%src(i, 1)
        iy = params%src(i, 2)
        do iz=1, params%nz
            write(file_unit, '(3I4, 9E11.3)') ix, iy, iz,                                                          &
                                              dx_source(ix, iy, iz), dy_source(ix, iy, iz), dz_source(ix, iy, iz), &
                                              field%hx(ix, iy, iz),  field%hy(ix, iy, iz),  field%hz(ix, iy, iz),  &
                                              ex_source(ix, iy, iz), ey_source(ix, iy, iz), ez_source(ix, iy, iz)
        end do
    end do
    
    close(file_unit)
end subroutine

end module
