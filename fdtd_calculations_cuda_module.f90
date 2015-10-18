module fdtd_calculations_cuda_module
use fdtd_data_cuda_module

implicit none

contains

attributes(global) subroutine update_h_field_cuda(params, field, run_num)
    !Input
    type(fdtd_params_cuda), intent(in) :: params
    type(fdtd_field_cuda),  intent(in) :: field
    integer, value,    intent(in) :: run_num

    !Local vars 
    integer                :: ix, iy, iz
    real, dimension(:,:,:) :: ex_source
    real, dimension(:,:,:) :: ey_source
    real, dimension(:,:,:) :: ez_source
    
    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z
    
    !Setup based on run_num 1..3
    if(run_num .eq. 1) then
        ex_source => field%ex3
        ey_source => field%ey3
        ez_source => field%ez3
    else if(run_num .eq. 2) then
        ex_source => field%ex1
        ey_source => field%ey1
        ez_source => field%ez1
    else
        ex_source => field%ex2
        ey_source => field%ey2
        ez_source => field%ez2
    end if
    
    !Update Hx
    if(ix >= 2 .and ix <= params%nx-1 .and.
       iy >= 1 .and iy <= params%ny-1 .and.
       iz >= 1 .and. iz<=params%nz-1) then
        field%hx(ix, iy, iz) = field%hx(ix, iy, iz) -                              &
                               params%dt/(params%mu_0 * params%dy) *               &
                               (ez_source(ix, iy+1, iz) - ez_source(ix, iy, iz)) + &
                               params%dt/(params%mu_0 * params%dz) *               &
                               (ey_source(ix, iy, iz+1) - ey_source(ix, iy, iz))
    end if
    
    !Update Hy
    if(ix >= 1 .and ix <= params%nx-1 .and.
       iy >= 2 .and iy <= params%ny-1 .and.
       iz >= 1 .and. iz<=params%nz-1) then
        field%hy(ix, iy, iz) = field%hy(ix, iy, iz) -                              &
                               params%dt/(params%mu_0 * params%dz) *               &
                               (ex_source(ix, iy, iz+1) - ex_source(ix, iy, iz)) + &
                               params%dt/(params%mu_0 * params%dx) *               &
                               (ez_source(ix+1, iy, iz) - ez_source(ix, iy, iz))
    end if
    
    !Update Hz
    if(ix >= 1 .and ix <= params%nx-1 .and.
       iy >= 1 .and iy <= params%ny-1 .and.
       iz >= 2 .and. iz<=params%nz-1) then
        field%hz(ix, iy, iz) = field%hz(ix, iy, iz) -                              &
                               params%dt/(params%mu_0 * params%dx) *               &
                               (ey_source(ix+1, iy, iz) - ey_source(ix, iy, iz)) + &
                               params%dt/(params%mu_0 * params%dy) *               &
                               (ex_source(ix, iy+1, iz) - ex_source(ix, iy, iz))
    end if
end subroutine


attributes(global) subroutine update_d_field_cuda(params, field, run_num)
    !Input
    type(fdtd_params_cuda), intent(in) :: params
    type(fdtd_field_cuda),  intent(in) :: field
    integer, value,    intent(in) :: run_num

    !Local vars
    integer                :: ix, iy, iz
    real, dimension(:,:,:) :: dx_target
    real, dimension(:,:,:) :: dy_target
    real, dimension(:,:,:) :: dz_target

    real, dimension(:,:,:) :: dx_source
    real, dimension(:,:,:) :: dy_source
    real, dimension(:,:,:) :: dz_source
    
    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    if(run_num .eq. 1) then
        dx_target => field%dx1
        dy_target => field%dy1
        dz_target => field%dz1
        
        dx_source => field%dx3
        dy_source => field%dy3
        dz_source => field%dz3
    else if(run_num .eq. 2) then
        dx_target => field%dx2
        dy_target => field%dy2
        dz_target => field%dz2
        
        dx_source => field%dx1
        dy_source => field%dy1        
        dz_source => field%dz1
    else
        dx_target => field%dx3
        dy_target => field%dy3
        dz_target => field%dz3

        dx_source => field%dx2        
        dy_source => field%dy2        
        dz_source => field%dz2
    end if
    
    !Update Dx
    if(ix >= 1 .and ix <= params%nx-1 .and.
       iy >= 2 .and iy <= params%ny-1 .and.
       iz >= 2 .and. iz<=params%nz-1) then
        dx_target(ix, iy, iz) = dx_source(ix, iy, iz) +                                                 &
                                params%dt/params%dy * (field%hz(ix, iy, iz) - field%hz(ix, iy-1, iz)) - &
                                params%dt/params%dz * (field%hy(ix, iy, iz) - field%hy(ix, iy, iz-1))
    end if
    
    !Update Dy
    if(ix >= 2 .and ix <= params%nx-1 .and.
       iy >= 1 .and iy <= params%ny-1 .and.
       iz >= 2 .and. iz<=params%nz-1) then
        dy_target(ix, iy, iz) = dy_source(ix, iy, iz) +                                                 &
                                params%dt/params%dz * (field%hx(ix, iy, iz) - field%hx(ix, iy, iz-1)) - &
                                params%dt/params%dx * (field%hz(ix, iy, iz) - field%hz(ix-1, iy, iz))
    end if
    
    !Update Dz
    if(ix >= 2 .and ix <= params%nx-1 .and.
       iy >= 2 .and iy <= params%ny-1 .and.
       iz >= 1 .and. iz<=params%nz-1) then
            dz_target(ix, iy, iz) = dz_source(ix, iy, iz) +                                                 &
                                    params%dt/params%dx * (field%hy(ix, iy, iz) - field%hy(ix-1, iy, iz)) - &
                                    params%dt/params%dy * (field%hx(ix, iy, iz) - field%hx(ix, iy-1, iz))    
    end if
end subroutine


attributes(global) subroutine update_e_field_cuda(params, field, run_num)
    !Input
    type(fdtd_params_cuda), intent(in) :: params
    type(fdtd_field_cuda),  intent(in) :: field
    integer, value     intent(in) :: run_num

    !Local vars
    integer                :: ix, iy, iz
    real, dimension(:,:,:) :: ex_target
    real, dimension(:,:,:) :: ey_target
    real, dimension(:,:,:) :: ez_target

    real, dimension(:,:,:) :: ex_source_1
    real, dimension(:,:,:) :: ex_source_2
    real, dimension(:,:,:) :: dx_source_1
    real, dimension(:,:,:) :: dx_source_2
    real, dimension(:,:,:) :: dx_source_3

    real, dimension(:,:,:) :: ey_source_1
    real, dimension(:,:,:) :: ey_source_2
    real, dimension(:,:,:) :: dy_source_1
    real, dimension(:,:,:) :: dy_source_2
    real, dimension(:,:,:) :: dy_source_3

    real, dimension(:,:,:) :: ez_source_1
    real, dimension(:,:,:) :: ez_source_2
    real, dimension(:,:,:) :: dz_source_1
    real, dimension(:,:,:) :: dz_source_2
    real, dimension(:,:,:) :: dz_source_3

    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    if(run_num .eq. 1) then
        ex_target => field%ex1
        ey_target => field%ey1
        ez_target => field%ez1
        
        ex_source_1 => field%ex3
        ex_source_2 => field%ex2
        dx_source_1 => field%dx1
        dx_source_2 => field%dx3
        dx_source_3 => field%dx2

        ey_source_1 => field%ey3
        ey_source_2 => field%ey2
        dy_source_1 => field%dy1
        dy_source_2 => field%dy3
        dy_source_3 => field%dy2

        ez_source_1 => field%ez3
        ez_source_2 => field%ez2
        dz_source_1 => field%dz1
        dz_source_2 => field%dz3
        dz_source_3 => field%dz2
    else if(run_num .eq. 2) then
        ex_target => field%ex2
        ey_target => field%ey2
        ez_target => field%ez2
        
        ex_source_1 => field%ex1
        ex_source_2 => field%ex3
        dx_source_1 => field%dx2
        dx_source_2 => field%dx1
        dx_source_3 => field%dx3

        ey_source_1 => field%ey1
        ey_source_2 => field%ey3
        dy_source_1 => field%dy2
        dy_source_2 => field%dy1
        dy_source_3 => field%dy3

        ez_source_1 => field%ez1
        ez_source_2 => field%ez3
        dz_source_1 => field%dz2
        dz_source_2 => field%dz1
        dz_source_3 => field%dz3
    else
        ex_target => field%ex3
        ey_target => field%ey3
        ez_target => field%ez3
        
        ex_source_1 => field%ex2
        ex_source_2 => field%ex1
        dx_source_1 => field%dx3
        dx_source_2 => field%dx2
        dx_source_3 => field%dx1

        ey_source_1 => field%ey2
        ey_source_2 => field%ey1
        dy_source_1 => field%dy3
        dy_source_2 => field%dy2
        dy_source_3 => field%dy1

        ez_source_1 => field%ez2
        ez_source_2 => field%ez1
        dz_source_1 => field%dz3
        dz_source_2 => field%dz2
        dz_source_3 => field%dz1
    end if
    
    !Update Ex
    if(ix >= 1 .and. ix <= params%nx-1 .and.
       iy >= 2 .and. iy <= params%ny-1 .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ex_target(ix, iy, iz) = (                                                                         &
                                 1/(2 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) + &
                                 2 * params%dt *                                                           &
                                 (                                                                        &  
                                  params%eps_0 * field%eps_s(ix, iy, iz) +                                 &
                                  field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                       &
                                 ) +                                                                      &
                                 field%sigma(ix, iy, iz) * params%dt * params%dt)                           &
                                ) *                                                                       &
                                (                                                                         &
                                 (                                                                        &
                                  4 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) +   &
                                  2 * params%dt *                                                          &
                                  (                                                                       &
                                   params%eps_0 * field%eps_s(ix, iy, iz) +                                &
                                   field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                      &
                                  ) -                                                                     &
                                  field%sigma(ix, iy, iz) * params%dt * params%dt                           &
                                 ) *                                                                      &
                                 ex_source_1(ix, iy, iz) -                                                &
                                 (2 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz)) *  & 
                                 ex_source_2(ix, iy, iz) +                                                &
                                 (2 * (params%dt + field%tau_d(ix, iy, iz))) * dx_source_1(ix, iy, iz) -   &
                                 (2 * params%dt + 4 * field%tau_d(ix, iy, iz)) * dx_source_2(ix, iy, iz) + &
                                 (2*field%tau_d(ix, iy, iz)) * dx_source_3(ix, iy, iz)                    &
                                )
    end if
    
    !Update Ey
    if(ix >= 2 .and. ix <= params%nx-1 .and.
       iy >= 1 .and. iy <= params%ny-1 .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ey_target(ix, iy, iz) = (                                                                         &
                                 1/(2 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) + &
                                 2 * params%dt *                                                           &
                                 (                                                                        &
                                  params%eps_0 * field%eps_s(ix, iy, iz) +                                 &
                                  field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                       &
                                 ) +                                                                      &
                                 field%sigma(ix, iy, iz) * params%dt * params%dt)                           &
                                ) *                                                                       &
                                (                                                                         &
                                 (                                                                        &
                                  4 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) +   &
                                  2 * params%dt *                                                          &  
                                  (                                                                       &
                                   params%eps_0 * field%eps_s(ix, iy, iz) +                                &
                                   field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                      &
                                  ) -                                                                     &
                                  field%sigma(ix, iy, iz) * params%dt * params%dt                           &
                                 ) *                                                                      &
                                 ey_source_1(ix, iy, iz) -                                                & 
                                 (2 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz)) *  &
                                 ey_source_2(ix, iy, iz) +                                                &
                                 (2 * (params%dt + field%tau_d(ix, iy, iz))) * dy_source_1(ix, iy, iz) -   &
                                 (2 * params%dt + 4 * field%tau_d(ix, iy, iz)) * dy_source_2(ix, iy, iz) + &
                                 (2 * field%tau_d(ix, iy, iz)) * dy_source_3(ix, iy, iz)                  &
                                )
    end if
    
    !Update Ez
    if(ix >= 2 .and. ix <= params%nx-1 .and.
       iy >= 1 .and. iy <= params%ny-1 .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ez_target(ix, iy, iz) = (                                                                         &
                                 1/(2 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) + &
                                 2 * params%dt *                                                           &
                                 (                                                                        &
                                  params%eps_0 * field%eps_s(ix, iy, iz) +                                 &
                                  field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                       &
                                 ) +                                                                      &
                                 field%sigma(ix, iy, iz) * params%dt * params%dt)                           &
                                ) *                                                                       &
                                (                                                                         &
                                 (                                                                        &
                                  4 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) +   &
                                  2 * params%dt *                                                          &  
                                  (                                                                       &
                                   params%eps_0 * field%eps_s(ix, iy, iz) +                                &
                                   field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                      &
                                  ) -                                                                     &
                                  field%sigma(ix, iy, iz) * params%dt * params%dt                           &
                                 ) *                                                                      &
                                 ez_source_1(ix, iy, iz) -                                                & 
                                 (2 * params%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz)) *  &
                                 ez_source_2(ix, iy, iz) +                                                &
                                 (2 * (params%dt + field%tau_d(ix, iy, iz))) * dz_source_1(ix, iy, iz) -   &
                                 (2 * params%dt + 4 * field%tau_d(ix, iy, iz)) * dz_source_2(ix, iy, iz) + &
                                 (2 * field%tau_d(ix, iy, iz)) * dz_source_3(ix, iy, iz)                  &
                                )
    end if
end subroutine


attributes(global) subroutine update_source_cuda(params, field, run_num, runs_count)
    !Input
    type(fdtd_params_cuda), intent(in) :: params
    type(fdtd_field_cuda),  intent(in) :: field
    integer, value,    intent(in) :: run_num
    integer, value,    intent(in) :: runs_count

    !Local vars 
    integer                :: ix, iy, iz
    integer                :: i
    real, dimension(:,:,:) :: dz_target
    real, dimension(:,:,:) :: dz_source
    integer                :: x, y, z

    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    if(run_num .eq. 1) then
        dz_target => field%dz1
        dz_source => field%dz3
    else if(run_num .eq. 2) then
        dz_target => field%dz2
        dz_source => field%dz1
    else
        dz_target => field%dz3
        dz_source => field%dz2
    end if

    !Update source
    if(ix == 0 .and. iy == 0 .and. iz == 0) then
        do i=1, params%nsrc
            x = params%src(i, 1)
            y = params%src(i, 2)
            z = params%src(i, 3)
    
            dz_target(x, y, z) = dz_source(x, y, z) +                                              &
                                 params%dt/params%dx * (field%hy(x, y, z) - field%hy(x-1, y, z)) - &
                                 params%dt/params%dy * (field%hx(x, y, z) - field%hx(x, y-1, z)) - &
                                 params%jz(((runs_count-1)*3)+1)
        end do
    end if
end subroutine


attributes(global) subroutine update_mur_boundary_cuda(params, field, run_num)
    !Input
    type(fdtd_params_cuda), intent(in) :: params
    type(fdtd_field_cuda),  intent(in) :: field
    integer, value,    intent(in) :: run_num
    
    !Local vars
    integer                :: ix, iy, iz
    real, dimension(:,:,:) :: ex_target
    real, dimension(:,:,:) :: ey_target
    real, dimension(:,:,:) :: ez_target

    real, dimension(:,:,:) :: ex_source
    real, dimension(:,:,:) :: ey_source
    real, dimension(:,:,:) :: ez_source
    
    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z
    
    !Setup based on run_num 1..3
    if(run_num .eq. 1) then
        ex_target => field%ex1
        ey_target => field%ey1
        ez_target => field%ez1
        
        ex_source => field%ex3
        ey_source => field%ey3
        ez_source => field%ez3
    else if(run_num .eq. 2) then
        ex_target => field%ex2
        ey_target => field%ey2
        ez_target => field%ez2
        
        ex_source => field%ex1
        ey_source => field%ey1
        ez_source => field%ez1
    else
        ex_target => field%ex3
        ey_target => field%ey3
        ez_target => field%ez3
        
        ex_source => field%ex2
        ey_source => field%ey2
        ez_source => field%ez2
    end if

    !Update Ex
    if(ix >= 1 .and. ix <= params%nx-1 .and.
       iy == 1 .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ex_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy, iz))) *  &
                                (                                                                                      &
                                 (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                 ex_target(ix, iy+1, iz) +                                                             &
                                 (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                 ex_source(ix, iy+1, iz) -                                                             &
                                 (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy,iz))) *    &
                                 ex_source(ix, iy, iz)                                                                 &
                                )
    end if

    if(ix >= 1 .and. ix <= params%nx-1 .and.
       iy == params%ny .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ex_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *  &
                                (                                                                                        &
                                 (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                 ex_target(ix, iy-1, iz) +                                                               &
                                 (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                 ex_source(ix, iy-1, iz) -                                                               &
                                 (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                                                   &
                                )
    end if

    if(ix >= 1 .and. ix <= params%nx-1 .and.
       iy >= 2 .and. iy <= params%ny-1 .and.
       iz == 1) then
        ex_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *  &
                                (                                                                                      &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                 ex_target(ix, iy, iz+1) +                                                             &
                                 (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                 ex_source(ix, iy, iz+1) -                                                             &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                                                 &
                                )
    end if

    if(ix >= 1 .and. ix <= params%nx-1 .and.
       iy >= 2 .and. iy <= params%ny-1 .and.
       iz == params%nz) then
        ex_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *  &
                                (                                                                                        &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                 ex_target(ix, iy, iz-1) +                                                               &
                                 (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                 ex_source(ix, iy, iz-1) -                                                               &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *   & 
                                 ex_source(ix, iy, iz)                                                                   &
                                )
    end if

    !Update Ey
    if(ix == 1 .and.
       iy >= 1 .and. iy <= params%ny-1 .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ey_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz))) *  &
                                (                                                                                      &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * &
                                 ey_target(ix+1, iy, iz) +                                                             &
                                 (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * &
                                 ey_source(ix+1, iy, iz) -                                                             &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                                                 &
                                )
    end if

    if(ix == params%nx .and.
       iy >= 1 .and. iy <= params%ny-1 .and.
       iz >= 2 .and. iz <= params%nz-1) then
        ey_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz)))  * &
                                (                                                                                        &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                 ey_source(ix-1, iy, iz) +                                                               &
                                 (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                 ey_source(ix-1, iy, iz) -                                                               &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                                                   &
                                )
    end if

    if(ix => 2 .and. ix <= params%nx-1 .and.
       iy >= 1 .and. iy <= params%ny-1 .and.
       iz == 1) then
        ey_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *  &
                                (                                                                                      &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                 ey_target(ix, iy,iz+1) +                                                              &
                                 (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                 ey_source(ix, iy, iz+1) -                                                             &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                                                 &
                                )
    end if

    if(ix => 2 .and. ix <= params%nx-1 .and.
       iy >= 1 .and. iy <= params%ny-1 .and.
       iz == params%nz) then
        ey_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *  &
                                (                                                                                        &
                                 (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                 ey_target(ix, iy, iz-1) +                                                               &
                                 (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                 ey_source(ix, iy, iz-1) -                                                               &
                                 (params%dt - params%dz *sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *    &     
                                 ey_source(ix, iy, iz)                                                                   &
                                )
    end if

    !Update Ez
    if(ix == 1 .and.
       iy >= 2 .and. iy <= params%ny-1 .and.
       iz => 1 .and iz <= params%nz-1) then
        ez_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz))) *  &
                                (                                                                                      &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * & 
                                 ez_target(ix+1, iy, iz) +                                                             &
                                 (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * &
                                 ez_source(ix+1, iy, iz) -                                                             & 
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz)))  *  &
                                 ez_source(ix, iy, iz)                                                                 &
                                )
    end if
      
    if(ix == params%nx .and.
       iy >= 2 .and. iy <= params%ny-1 .and.
       iz => 1 .and iz <= params%nz-1) then
        ez_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz))) *  &
                                (                                                                                        &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                 ez_target(ix-1, iy, iz) +                                                               &
                                 (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                 ez_source(ix-1, iy, iz) -                                                               &
                                 (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz))) *   &
                                 ez_source(ix, iy, iz)                                                                   &
                                )
    end if
    
    if(ix >= 2 .and. ix <= params%nx-1 .and.
       iy == 1 .and.
       iz => 1 .and iz <= params%nz-1) then 
        ez_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy, iz))) *  &
                                (                                                                                      &
                                 (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                 ez_target(ix, iy+1, iz) +                                                             & 
                                 (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                 ez_source(ix, iy+1, iz) -                                                             &
                                 (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy, iz))) *   &
                                 ez_source(ix, iy, iz)                                                                 &
                                )
    end if
      
    if(ix >= 2 .and. ix <= params%nx-1 .and.
       iy == params%ny .and.
       iz => 1 .and iz <= params%nz-1) then 
            ez_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *  &
                                    (                                                                                        &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                     ez_target(ix, iy-1, iz) +                                                               &
                                     (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                     ez_source(ix, iy-1, iz) -                                                               &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *   &
                                     ez_source(ix, iy, iz)                                                                   &
                                    )
    end if
end

end module