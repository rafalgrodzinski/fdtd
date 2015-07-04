module fdtd_data_module
use utils_module

implicit none

type fdtd_state
    !read from file
    integer :: nx, ny, nz
    integer :: t_max
    integer :: nf
    character(len=128) :: rdirname
    character(len=32) :: rprefix
    character(len=128) :: dirname
    character :: outformat*3, chr
    integer :: impulse_resp_flag, pec_flag, read_env_flag
    integer :: output_flag, bzip2_flag
    integer, dimension(0:2) :: output_start, output_finish
    character :: stype*4
    integer :: elem_lambda
    real :: w_freq
    real :: pwidth, pmodufreq
    integer :: nsrc
    integer, allocatable :: src(:,:)
    integer :: ptype
    real :: fsigma, feps_s, feps_i, ftau_d
    real :: feps_m, ftau_2
    real :: a
    
    !generated
    real :: pi = acos(-1.0) !Delicious pie
    real :: c = 3.0*10.0**8 !Light speed (v in m/s)
    real :: timeskip !Time step skip
    real :: lambda !Wave length (meters)
    real :: dt !Length of the time step
    real :: dx, dy, dz !Distance between 2 cells
    real :: mu_0 !mu0, permeability of free space (in henry/meter)
    real :: eps_0 !Epsilon0, permittivity of free space (in farad/meter)
    
    real, dimension(:), allocatable :: jz
end type


type fdtd_field
    real, dimension(:,:,:), pointer :: ex1,ex2,ex3
    real, dimension(:,:,:), pointer :: ey1,ey2,ey3 
    real, dimension(:,:,:), pointer :: ez1,ez2,ez3 

    real, dimension(:,:,:), pointer :: hx
    real, dimension(:,:,:), pointer :: hy
    real, dimension(:,:,:), pointer :: hz

    real, dimension(:,:,:), pointer :: dx1,dx2,dx3 
    real, dimension(:,:,:), pointer :: dy1,dy2,dy3 
    real, dimension(:,:,:), pointer :: dz1,dz2,dz3 

    real, dimension(:,:,:), pointer :: eps_i, eps_s
    real, dimension(:,:,:), pointer :: tau_d, sigma

    !mur boundary
    real, dimension(:,:,:), pointer :: rp_x_1, rp_x_end
    real, dimension(:,:,:), pointer :: rp_y_1, rp_y_end
    real, dimension(:,:,:), pointer :: rp_z_1, rp_z_end
end type


contains

subroutine init_fdtd_state(state, file_name)
    !input
    type(fdtd_state), pointer, intent(inout) :: state
    character(len=*), intent(in)             :: file_name
    !local vars
    integer, parameter :: file_unit = 51
    integer            :: error_code
    character          :: dummy_char
    integer            :: i
    
    allocate(state)
    
    open(unit=file_unit, file=file_name, status="old", iostat=error_code)
    call check_error(error_code, "Couldn't open file "//file_name)

    read(51, *, iostat=error_code) dummy_char, state%nx, state%ny, state%nz
    read(51, *, iostat=error_code) dummy_char, state%t_max
    read(51, *, iostat=error_code) dummy_char, state%nf
    read(51, *, iostat=error_code) dummy_char, state%rdirname
    read(51, *, iostat=error_code) dummy_char, state%rprefix
    read(51, *, iostat=error_code) dummy_char, state%dirname
    read(51, *, iostat=error_code) dummy_char, state%outformat
    read(51, *, iostat=error_code) dummy_char, state%impulse_resp_flag
    read(51, *, iostat=error_code) dummy_char, state%pec_flag
    read(51, *, iostat=error_code) dummy_char, state%read_env_flag
    read(51, *, iostat=error_code) dummy_char, state%output_flag
    read(51, *, iostat=error_code) dummy_char, state%bzip2_flag
    read(51, *, iostat=error_code) dummy_char, state%output_start(0:2)
    read(51, *, iostat=error_code) dummy_char, state%output_finish(0:2)
    read(51, *, iostat=error_code) dummy_char, state%stype
    read(51, *, iostat=error_code) dummy_char, state%elem_lambda
    read(51, *, iostat=error_code) dummy_char, state%w_freq
    read(51, *, iostat=error_code) dummy_char, state%pwidth
    read(51, *, iostat=error_code) dummy_char, state%pmodufreq
    read(51, *, iostat=error_code) dummy_char, state%nsrc

    allocate (state%src(1:state%nsrc, 1:3))

    do i = 1, state%nsrc
      read(51, *, iostat=error_code) dummy_char, state%src(i,1:3)
    enddo

    read(51, *, iostat=error_code) dummy_char, state%ptype

    read(51, *, iostat=error_code) dummy_char, state%fsigma
    read(51, *, iostat=error_code) dummy_char, state%feps_s
    read(51, *, iostat=error_code) dummy_char, state%feps_i
    read(51, *, iostat=error_code) dummy_char, state%ftau_d
    read(51, *, iostat=error_code) dummy_char, state%feps_m
    read(51, *, iostat=error_code) dummy_char, state%ftau_2

    read(51, *, iostat=error_code) dummy_char, state%A
    
    close(file_unit)
    
    !generate rest of the values
    state%timeskip = 1.0
    state%lambda = state%c / state%w_freq
    state%dx = state%lambda/state%elem_lambda
    state%dy = state%dx
    state%dz = state%dx
    state%dt = 1.0d0 * state%timeskip / (state%c * sqrt(1.0d0/(state%dx**2) + 1.0d0/(state%dy**2) + 1.0d0/(state%dz**2)))
    state%mu_0 = 4*state%pi*10**(-7.0)
    state%eps_0 = 1.0/(state%mu_0 * state%c * state%c)
end


subroutine delete_fdtd_state(state)
    type(fdtd_state), pointer, intent(inout) :: state
    
    deallocate(state%src)
    deallocate(state%jz)
    deallocate(state)
end


subroutine init_fdtd_field(field, state)
    !input
    type(fdtd_field), pointer, intent(inout) :: field
    type(fdtd_state), pointer, intent(in)    :: state
    !local vars
    integer :: ix, iy, iz
    
    allocate(field)
    
    allocate(field%ex1(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%ex2(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%ex3(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%ey1(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%ey2(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%ey3(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%ez1(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%ez2(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%ez3(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%hx(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%hy(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%hz(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%dx1(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%dx2(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%dx3(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%dy1(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%dy2(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%dy3(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%dz1(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%dz2(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%dz3(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%eps_i(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%eps_s(1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%tau_d(1:state%nx, 1:state%ny, 1:state%nz))
    allocate(field%sigma(1:state%nx, 1:state%ny, 1:state%nz))
    
    field%ex1(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%ex2(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%ex3(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%ey1(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%ey2(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%ey3(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%ez1(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%ez2(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%ez3(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%hx(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%hy(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%hz(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%dx1(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%dx2(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%dx3(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%dy1(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%dy2(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%dy3(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%dz1(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%dz2(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%dz3(1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    field%eps_i(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%eps_s(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%tau_d(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    field%sigma(1:state%nx, 1:state%ny, 1:state%nz) = 0.0
    
    field%sigma(1:state%nx, 1:state%ny, 1:state%nz) = state%fsigma
    field%eps_s(1:state%nx, 1:state%ny, 1:state%nz) = state%feps_s
    field%eps_i(1:state%nx, 1:state%ny, 1:state%nz) = state%feps_i
    field%tau_d(1:state%nx, 1:state%ny, 1:state%nz) = state%ftau_d
    
    !Initialise mur boundary
    allocate(field%rp_x_1(1:2, 1:state%ny, 1:state%nz))
    allocate(field%rp_x_end(state%nx-1:state%nx, 1:state%ny, 1:state%nz))

    allocate(field%rp_y_1(1:state%nx, 1:2, 1:state%nz))
    allocate(field%rp_y_end(1:state%nx, state%ny-1:state%ny, 1:state%nz))

    allocate(field%rp_z_1(1:state%nx, 1:state%ny, 1:2))
    allocate(field%rp_z_end(1:state%nx, 1:state%ny, state%nz-1:state%nz))
    
    !Setup rp_x
    field%rp_x_1(1:2, 1:state%ny, 1:state%nz) = 0.0
    field%rp_x_end(state%nx-1:state%nx, 1:state%ny, 1:state%nz) = 0.0

    do iz=1, state%nz
        do iy=1, state%ny
            do ix=1, 2
                field%rp_x_1(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                    &
                                            (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                            &
                                            (1 + cmplx(0.0, 2 * state%pi * state%w_freq * field%tau_d(ix, iy,iz))) -         &
                                            cmplx(0, field%sigma(ix, iy, iz) / (2 * state%pi * state%w_freq * state%eps_0)))
            end do
	    end do
    end do

    do iz=1, state%nz
	    do iy=1, state%ny
            do ix=state%nx-1, state%nx 
                field%rp_x_end(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                               &
                                                  (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                   &
                                                  (1 + cmplx(0, 2 * state%pi * state%w_freq * field%tau_d(ix, iy, iz))) - &
                                                  cmplx(0, field%sigma(ix, iy, iz) / (2 * state%pi * state%w_freq * state%eps_0)))
            end do
	    end do
    end do
    
    !Setup rp_y
    field%rp_y_1(1:state%nx, 1:2, 1:state%nz) = 0.0
    field%rp_y_end(1:state%nx, state%ny-1:state%ny, 1:state%nz) = 0.0

    do iz=1, state%nz
	    do iy=1, 2 
            do ix=1, state%nx 
                field%rp_y_1(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                                &
                                                (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                    &
                                                (1 + cmplx(0, 2 * state%pi * state%w_freq * field%tau_d(ix, iy, iz))) -  &
                                                cmplx(0, field%sigma(ix, iy, iz) / (2 * state%pi * state%w_freq * state%eps_0)))
            end do
	    end do
    end do


    do iz=1, state%nz
	    do iy=state%ny-1, state%ny
            do ix=1, state%nx 
                field%rp_y_end(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                               &
                                                  (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                   &
                                                  (1 + cmplx(0, 2 * state%pi * state%w_freq * field%tau_d(ix, iy, iz))) - &
                                                  cmplx(0, field%sigma(ix, iy, iz) / (2 * state%pi * state%w_freq * state%eps_0)))
            end do
	    end do
    end do
    
    !Setup rp_z
    field%rp_z_1(1:state%nx, 1:state%ny, 1:2) = 0.0
    field%rp_z_end(1:state%nx, 1:state%ny, state%nz-1:state%nz) = 0.0

    do iz=1, 2 
	    do iy=1, state%ny
            do ix=1, state%nx
                field%rp_z_1(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                               &
                                                (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                   &
                                                (1 + cmplx(0, 2 * state%pi * state%w_freq * field%tau_d(ix, iy, iz))) - &
                                                cmplx(0, field%sigma(ix, iy, iz) / (2 * state%pi * state%w_freq * state%eps_0)))
            end do
	    end do
    end do

    do iz=state%nz-1, state%nz 
    	do iy=1, state%ny
            do ix=1, state%nx
                field%rp_z_end(ix, iy, iz) = real(field%eps_i(ix, iy, iz) +                                               &
                                                  (field%eps_s(ix, iy, iz) - field%eps_i(ix, iy, iz)) /                   &
                                                  (1 + cmplx(0, 2 * state%pi * state%w_freq * field%tau_d(ix, iy, iz))) - &
                                                  cmplx(0, field%sigma(ix, iy, iz) / (2 * state%pi * state%w_freq * state%eps_0)))
            end do
	    end do
    end do
end


subroutine delete_fdtd_field(field)
    !input
    type(fdtd_field), pointer, intent(inout) :: field
    
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
end


subroutine load_materials(state, field, specs_file_name, materials_path)
    !input
    type(fdtd_state), pointer, intent(in) :: state
    type(fdtd_field), pointer, intent(in) :: field  
    character(len=*), intent(in)          :: specs_file_name
    character(len=*), intent(in)          :: materials_path
    !local vars
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
    do iz=1, state%nz
        !generte file name, starting with v1_00001.pgm
        write(material_file_name, fmt='(I5)'), iz
        material_file_name = "0000" // adjustl(material_file_name)
        material_file_name = material_file_name(len(trim(material_file_name))-4 : len(trim(material_file_name)))
        material_file_name = materials_path // "v1_" // trim(material_file_name) // ".pgm"
    
        open(unit=file_unit, file=material_file_name, status="old", iostat=error_code)
        call check_error(error_code, "Couldn't open file "//material_file_name)
        
        read(file_unit, *) dummy_char, dummy_char, dummy_char, material_width, material_height, dummy_char
        
        do iy=1, state%ny
            do ix=1, state%nx
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
end


subroutine setup_source(state, field)
    !input
    type(fdtd_state), pointer, intent(in) :: state
    type(fdtd_field), pointer, intent(in) :: field
    !local vars
    integer :: fine
    integer :: temp
    integer :: i, i2
    integer :: istart
    real, dimension(:), allocatable :: tmpdata, tmpdata2
    
    allocate(state%jz(1:2**16))
    allocate(tmpdata(-2**16:2**16))
    allocate(tmpdata2(-2**16:2**16))
    
    !fine & temp
    fine = int(2**13 * state%pwidth * state%w_freq * state%dt)

    temp = 1/(state%pwidth * state%w_freq)/(state%dt / fine)/2
    
    !tmpdata
    do i=-2**14, 2**14
        tmpdata(i)=cmplx(0.0)
    end do
      
    do i=-temp-1,temp+1
         tmpdata(i)= exp(-(real(i)/((real(temp)+1.0)/4.0))**2)
    end do
      
    do i = -temp-1,temp+1
      tmpdata(i) = tmpdata(i) *cos(2.0*acos(-1.0)* state%pmodufreq * state%w_freq*i*(state%dt/fine))     
    enddo

    !istart
    do i=-2**12,2**12-1
         if( (abs(tmpdata(i)).gt.1e-9).and.(mod(i,fine).eq.0)) then
            istart = i 
            exit
         endif
    enddo
    
    !setup jz 1/2
    state%jz = 0.0
    
    i2 = 0
    do i=istart, temp+1, fine
          i2=i2+1
          state%jz(i2)=tmpdata(i)*10**(-15.)/state%dt/3.0
    end do
    
    !setup tmpdata2
    do i=2, 2**14
        tmpdata2(i-1)=(real((state%jz(i+1)-state%jz(i))/state%dt)+real((state%jz(i)-state%jz(i-1))/state%dt)) &
        /2.0*(state%dt*state%dz)/(state%dx*state%dy*state%dz)
    end do
    
    !setup jz 2/2
    do i=1, 2**14
        state%jz(i) = tmpdata2(i)
    end do
end


subroutine write_result(state, field, run_num, runs_count, output_path)
    !input
    type(fdtd_state), pointer, intent(in) :: state
    type(fdtd_field), pointer, intent(in) :: field
    integer, intent(in)                   :: run_num
    integer, intent(in)                   :: runs_count
    character(len=*), intent(in)          :: output_path
    !local vars
    integer, parameter :: file_unit = 100
    character(len=128) :: output_file_name
    integer            :: error_code
    integer            :: i
    integer            :: ix, iy, iz, count
    real, dimension(:,:,:), pointer :: ex_source
    real, dimension(:,:,:), pointer :: ey_source
    real, dimension(:,:,:), pointer :: ez_source
    
    real, dimension(:,:,:), pointer :: dx_source
    real, dimension(:,:,:), pointer :: dy_source
    real, dimension(:,:,:), pointer :: dz_source
    
    !setup based on run_num 1..3
    if(run_num .eq. 1) then
        ex_source => field%ex1
        ey_source => field%ey1
        ez_source => field%ez1
        
        dx_source => field%dx1
        dy_source => field%dy1
        dz_source => field%dz1
    else if(run_num .eq. 2) then
        ex_source => field%ex2
        ey_source => field%ey2
        ez_source => field%ez2
        
        dx_source => field%dx2
        dy_source => field%dy2
        dz_source => field%dz2
    else
        ex_source => field%ex3
        ey_source => field%ey3
        ez_source => field%ez3
        
        dx_source => field%dx3
        dy_source => field%dy3
        dz_source => field%dz3
    end if
    
    !Output x
    !generte file name, starting with E_field_x_00001.out
    write(output_file_name, fmt='(I5)'), runs_count
    output_file_name = "0000" // adjustl(output_file_name)
    output_file_name = output_file_name(len(trim(output_file_name))-4 : len(trim(output_file_name)))
    output_file_name = output_path // "E_field_x_" // trim(output_file_name) // ".out"
    
    open(file_unit, file=output_file_name, status="new", access="sequential", form="formatted", iostat=error_code)
    call check_error(error_code, "Couldn't open file " // output_file_name)
    
    do i=1, state%nsrc
        iy = state%src(i, 2)
        iz = state%src(i, 3)
        do ix=1, state%nx
            write(file_unit, '(3I4, 9E11.3)') ix, iy, iz,                                                          &
                                              dx_source(ix, iy, iz), dy_source(ix, iy, iz), dz_source(ix, iy, iz), &
                                              field%hx(ix, iy, iz),  field%hy(ix, iy, iz),  field%hz(ix, iy, iz),  &
                                              ex_source(ix, iy, iz), ey_source(ix, iy, iz), ez_source(ix, iy, iz)
        end do
    end do
    
    close(file_unit)
    
    !Output y
    !generte file name, starting with E_field_y_00001.out
    write(output_file_name, fmt='(I5)'), runs_count
    output_file_name = "0000" // adjustl(output_file_name)
    output_file_name = output_file_name(len(trim(output_file_name))-4 : len(trim(output_file_name)))
    output_file_name = output_path // "E_field_y_" // trim(output_file_name) // ".out"
    
    open(file_unit, file=output_file_name, status="new", access="sequential", form="formatted", iostat=error_code)
    call check_error(error_code, "Couldn't open file " // output_file_name)
    
    do i=1, state%nsrc
        ix = state%src(i, 1)
        iz = state%src(i, 3)
        do iy=1, state%ny
            write(file_unit, '(3I4, 9E11.3)') ix, iy, iz,                                                          &
                                              dx_source(ix, iy, iz), dy_source(ix, iy, iz), dz_source(ix, iy, iz), &
                                              field%hx(ix, iy, iz),  field%hy(ix, iy, iz),  field%hz(ix, iy, iz),  &
                                              ex_source(ix, iy, iz), ey_source(ix, iy, iz), ez_source(ix, iy, iz)
        end do
    end do
    
    close(file_unit)
    
    !Output z
    !generte file name, starting with E_field_z_00001.out
    write(output_file_name, fmt='(I5)'), runs_count
    output_file_name = "0000" // adjustl(output_file_name)
    output_file_name = output_file_name(len(trim(output_file_name))-4 : len(trim(output_file_name)))
    output_file_name = output_path // "E_field_z_" // trim(output_file_name) // ".out"
    
    open(file_unit, file=output_file_name, status="new", access="sequential", form="formatted", iostat=error_code)
    call check_error(error_code, "Couldn't open file " // output_file_name)
    
    do i=1, state%nsrc
        ix = state%src(i, 1)
        iy = state%src(i, 2)
        do iz=1, state%nz
            write(file_unit, '(3I4, 9E11.3)') ix, iy, iz,                                                          &
                                              dx_source(ix, iy, iz), dy_source(ix, iy, iz), dz_source(ix, iy, iz), &
                                              field%hx(ix, iy, iz),  field%hy(ix, iy, iz),  field%hz(ix, iy, iz),  &
                                              ex_source(ix, iy, iz), ey_source(ix, iy, iz), ez_source(ix, iy, iz)
        end do
    end do
    
    close(file_unit)
end

end