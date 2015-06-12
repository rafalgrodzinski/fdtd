include "utils_module.f90"


module fdtd_data_module
use utils_module

implicit none

type fdtd_state
    integer :: nx,ny,nz
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
    real :: A
end type


type fdtd_field
    real, dimension(:,:,:), allocatable :: ex1,ex2,ex3
    real, dimension(:,:,:), allocatable :: ey1,ey2,ey3 
    real, dimension(:,:,:), allocatable :: ez1,ez2,ez3 

    real, dimension(:,:,:), allocatable :: hx
    real, dimension(:,:,:), allocatable :: hy
    real, dimension(:,:,:), allocatable :: hz

    real, dimension(:,:,:), allocatable :: dx1,dx2,dx3 
    real, dimension(:,:,:), allocatable :: dy1,dy2,dy3 
    real, dimension(:,:,:), allocatable :: dz1,dz2,dz3 

    real, dimension(:,:,:), allocatable :: eps_i, eps_s
    real, dimension(:,:,:), allocatable :: tau_d, sigma

    real, dimension(:,:,:), allocatable :: rp_x_1, rp_x_end
    real, dimension(:,:,:), allocatable :: rp_y_1, rp_y_end

    real, dimension(:,:,:), allocatable :: rp_z_1, rp_z_end
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

    allocate (state%src(0:state%nsrc-1, 0:2))

    do i = 0, state%nsrc-1
      read(51, *, iostat=error_code) dummy_char, state%src(i,0:2)
    enddo

    read(51, *, iostat=error_code) dummy_char, state%ptype

    read(51, *, iostat=error_code) dummy_char, state%fsigma
    read(51, *, iostat=error_code) dummy_char, state%feps_s
    read(51, *, iostat=error_code) dummy_char, state%feps_i
    read(51, *, iostat=error_code) dummy_char, state%ftau_d
    read(51, *, iostat=error_code) dummy_char, state%feps_m
    read(51, *, iostat=error_code) dummy_char, state%ftau_2

    read(51, *, iostat=error_code) dummy_char, state%A
end


subroutine delete_fdtd_state(state)
    type(fdtd_state), pointer, intent(inout) :: state
    
    deallocate(state)
end


subroutine init_fdtd_field(field, state)
    type(fdtd_field), pointer, intent(inout) :: field
    type(fdtd_state), pointer, intent(in)    :: state
end


subroutine delete_fdtd_field(field)
    type(fdtd_field), pointer, intent(inout) :: field
end

end