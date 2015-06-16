module fdtd_calculations_module
use fdtd_data_module
use utils_module

implicit none

contains

subroutine update_h_field(field, state, run_num)
    !input
    type(fdtd_field), pointer, intent(inout) :: field
    type(fdtd_state), pointer, intent(in)    :: state
    integer, intent(in)                      :: run_num
    !local vars
    integer :: ix, iy, iz
    real, dimension(:,:,:), pointer :: ex_field
    real, dimension(:,:,:), pointer :: ey_field
    real, dimension(:,:,:), pointer :: ez_field
    
    if(run_num .eq. 0) then
        ex_field => field%ex2
        ey_field => field%ey2
        ez_field => field%ez2
    else if(run_num .eq. 1) then
        ex_field => field%ex3
        ey_field => field%ey3
        ez_field => field%ez3
    else
        ex_field => field%ex1
        ey_field => field%ey1
        ez_field => field%ez1
    end if
    
    !Update Hx
    do iz=1, state%nz-1
        do iy=1, state%ny-1
	        do ix=2, state%nx-1
                field%hx(ix, iy, iz) = field%hx(ix, iy, iz) -             &
                                       state%dt/(state%mu_0 * state%dy) * &
                                       (ez_field(ix, iy+1, iz) -          &
                                       ez_field(ix, iy, iz)) +            &
                                       state%dt/(state%mu_0 * state%dz) * &
                                       (ey_field(ix, iy, iz+1) -          &
                                       ey_field(iz, iy, iz))
            end do
	    end do
    end do
    
    !Update Hy
    do iz=1, state%nz-1
        do iy=2, state%ny-1
	        do ix=1, state%nx-1
                field%hy(ix, iy, iz) = field%hy(ix, iy, iz) -             &
                                       state%dt/(state%mu_0 * state%dz) * &
                                       (ex_field(ix, iy, iz+1) -          &
                                       ex_field(ix, iy, iz)) +            &
                                       state%dt/(state%mu_0 * state%dx) * &
                                       (ez_field(ix+1, iy, iz) -          &
                                       ez_field(ix, iy, iz))
            end do
	    end do
    end do
    
    !Update Hz
    do iz=2, state%nz-1
        do iy=1, state%ny-1
	        do ix=1, state%nx-1
                field%hz(ix, iy, iz) = field%hz(ix, iy, iz) -             &
                                       state%dt/(state%mu_0 * state%dx) * &
                                       (ey_field(ix+1, iy, iz) -          &
                                       ey_field(ix, iy, iz)) +            &
                                       state%dt/(state%mu_0 * state%dy) * &
                                       (ex_field(ix, iy+1, iz) -          &
                                       ex_field(ix, iy, iz))
            end do
	    end do
    end do
end


subroutine update_d_field(field, state, run_num)
    !input
    type(fdtd_field), pointer, intent(inout) :: field
    type(fdtd_state), pointer, intent(in)    :: state
    integer, intent(in)                      :: run_num
    !local vars
    integer :: ix, iy, iz
    real, dimension(:,:,:), pointer :: dx_source
    real, dimension(:,:,:), pointer :: dx_target
    real, dimension(:,:,:), pointer :: dy_source
    real, dimension(:,:,:), pointer :: dy_target
    real, dimension(:,:,:), pointer :: dz_source
    real, dimension(:,:,:), pointer :: dz_target
    
    if(run_num .eq. 0) then
        dx_source => field%dx2
        dx_target => field%dx1
        dy_source => field%dx2
        dy_target => field%dx1
        dz_source => field%dx2
        dz_target => field%dx1
    else if(run_num .eq. 1) then
    else
    end if
    
    !Update Dx
    do iz=2, state%nz-1
        do iy=2, state%ny-1
	        do ix=1, state%nx-1
                dx_target(ix, iy, iz) = dx_source(ix, iy, iz) +                     &
                                        state%dt/state%dy * (field%hz(ix, iy, iz) - &
                                        field%hz(ix, iy-1, iz)) -                   &
                                        state%dt/state%dz * (field%hy(ix, iy, iz) - &
                                        field%hy(ix, iy, iz-1))
            end do
	    end do
    end do
    
    !Update Dy
    do iz=2, state%nz-1
        do iy=1, state%ny-1
	        do ix=2, state%nx-1
                dy_target(ix, iy, iz) = dy_source(ix, iy, iz) +                     &
                                        state%dt/state%dz * (field%hx(ix, iy, iz) - &
                                        field%hx(ix, iy, iz-1)) -                   &
                                        state%dt/state%dx * (field%hz(ix, iy, iz) - &
                                        field%hz(ix-1, iy, iz))
            end do
	    end do
    end do
    
    !Update Dz
    do iz=1, state%nz-1
        do iy=2, state%ny-1
            do ix=2, state%nx-1
                dz_target(ix, iy, iz) = dz_source(ix, iy, iz) +                     &
                                        state%dt/state%dx * (field%hy(ix, iy, iz) - &
                                        field%hy(ix-1, iy, iz)) -                   &
                                        state%dt/state%dy * (field%hx(ix, iy, iz) - &
                                        field%hx(ix, iy-1, iz))
            end do
        end do
    end do
end


subroutine update_e_field(field)
    type(fdtd_field), pointer, intent(inout) :: field
    
    !Update Ex
    
    !Update Ey
    
    !Update Ez
end

end