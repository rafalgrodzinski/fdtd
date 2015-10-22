module fdtd_calculations_cuda_module

use cudafor

use fdtd_data_cuda_module


implicit none

contains

attributes(global) subroutine update_h_field_cuda(hx, hy, hz,
                                                  ex_source, ey_source, ez_source,
                                                  nx, ny, nz,
                                                  dt, dx, dy, dz,
                                                  mu_0)

    !Input
    real, dimension(:,:,:), intent(in) :: hx, hy, hz
    real, dimension(:,:,:), intent(in) :: ex_source, ey_source, ez_source
    integer, intent(in)                :: nx, ny, nz
    real, intent(in)                   :: dt, dx, dy, dz
    real, intent(in)                   :: mu_0
    !type(fdtd_params_cuda), intent(in) :: params
    !type(fdtd_field_cuda), intent(in)  :: field
    !integer, value, intent(in)         :: run_num

    !Local vars 
    integer                         :: ix, iy, iz
    !real, dimension(:,:,:), pointer :: ex_source
    !real, dimension(:,:,:), pointer :: ey_source
    !real, dimension(:,:,:), pointer :: ez_source
    
    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z
    
    !Setup based on run_num 1..3
    !if(run_num .eq. 1) then
    !    ex_source => field%ex3
    !    ey_source => field%ey3
    !    ez_source => field%ez3
    !else if(run_num .eq. 2) then
    !    ex_source => field%ex1
    !    ey_source => field%ey1
    !    ez_source => field%ez1
    !else
    !    ex_source => field%ex2
    !    ey_source => field%ey2
    !    ez_source => field%ez2
    !end if
    
    !Update Hx
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 1 .and. iz<=nz-1) then
        hx(ix, iy, iz) = hx(ix, iy, iz) -                                          &
                               dt/(mu_0 * dy) *                                    &
                               (ez_source(ix, iy+1, iz) - ez_source(ix, iy, iz)) + &
                               dt/(mu_0 * dz) *                                    &
                               (ey_source(ix, iy, iz+1) - ey_source(ix, iy, iz))
    end if
    
    !Update Hy
    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz >= 1 .and. iz<=nz-1) then
        hy(ix, iy, iz) = hy(ix, iy, iz) -                                          &
                               dt/(mu_0 * dz) *                                    &
                               (ex_source(ix, iy, iz+1) - ex_source(ix, iy, iz)) + &
                               dt/(mu_0 * dx) *                                    &
                               (ez_source(ix+1, iy, iz) - ez_source(ix, iy, iz))
    end if
    
    !Update Hz
    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz<=nz-1) then
        hz(ix, iy, iz) = hz(ix, iy, iz) -                                          &
                               dt/(mu_0 * dx) *                                    &
                               (ey_source(ix+1, iy, iz) - ey_source(ix, iy, iz)) + &
                               dt/(mu_0 * dy) *                                    &
                               (ex_source(ix, iy+1, iz) - ex_source(ix, iy, iz))
    end if
end subroutine


attributes(global) subroutine update_d_field_cuda(dx_target, dy_target, dz_target,
                                                  dx_source, dy_source, dz_source,
                                                  hx, hy, hz,
                                                  nx, ny, nz,
                                                  dt, dx, dy, dz)

    !Input
    real, dimension(:,:,:), intent(in) :: dx_target, dy_target, dz_target
    real, dimension(:,:,:), intent(in) :: dx_source, dy_source, dz_source
    real, dimension(:,:,:), intent(in) :: hx, hy, hz
    integer, intent(in)                :: nx, ny, nz
    real, intent(in)                   :: dt, dx, dy, dz
    !type(fdtd_params_cuda), intent(in) :: params
    !type(fdtd_field_cuda), intent(in)  :: field
    !integer, value, intent(in)         :: run_num

    !Local vars
    integer                         :: ix, iy, iz
    !real, dimension(:,:,:), pointer :: dx_target
    !real, dimension(:,:,:), pointer :: dy_target
    !real, dimension(:,:,:), pointer :: dz_target

    !real, dimension(:,:,:), pointer :: dx_source
    !real, dimension(:,:,:), pointer :: dy_source
    !real, dimension(:,:,:), pointer :: dz_source
    
    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    !if(run_num .eq. 1) then
    !    dx_target => field%dx1
    !    dy_target => field%dy1
    !    dz_target => field%dz1
        
    !    dx_source => field%dx3
    !    dy_source => field%dy3
    !    dz_source => field%dz3
    !else if(run_num .eq. 2) then
    !    dx_target => field%dx2
    !    dy_target => field%dy2
    !    dz_target => field%dz2
        
    !    dx_source => field%dx1
    !    dy_source => field%dy1        
    !    dz_source => field%dz1
    !else
    !    dx_target => field%dx3
    !    dy_target => field%dy3
    !    dz_target => field%dz3

    !    dx_source => field%dx2        
    !    dy_source => field%dy2        
    !    dz_source => field%dz2
    !end if
    
    !Update Dx
    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz<=nz-1) then
        dx_target(ix, iy, iz) = dx_source(ix, iy, iz) +                       &
                                dt/dy * (hz(ix, iy, iz) - hz(ix, iy-1, iz)) - &
                                dt/dz * (hy(ix, iy, iz) - hy(ix, iy, iz-1))
    end if
    
    !Update Dy
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz<=nz-1) then
        dy_target(ix, iy, iz) = dy_source(ix, iy, iz) +                       &
                                dt/dz * (hx(ix, iy, iz) - hx(ix, iy, iz-1)) - &
                                dt/dx * (hz(ix, iy, iz) - hz(ix-1, iy, iz))
    end if
    
    !Update Dz
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz >= 1 .and. iz<=nz-1) then
            dz_target(ix, iy, iz) = dz_source(ix, iy, iz) +                       &
                                    dt/dx * (hy(ix, iy, iz) - hy(ix-1, iy, iz)) - &
                                    dt/dy * (hx(ix, iy, iz) - hx(ix, iy-1, iz))    
    end if
end subroutine


attributes(global) subroutine update_e_field_cuda(ex_target, ey_target, ez_target, &
                                                  ex_source_1, ey_source_1, ez_source_1, &
                                                  ex_source_2, ex_source_2, ex_source_2, &
                                                  dx_source_1, dy_source_1, dz_source_1, &
                                                  dx_source_2, dy_source_2, dz_source_2, &
                                                  dx_source_3, dy_source_3, dz_source_3, &
                                                  eps_i, eps_s, &
                                                  tau_d, sigma, &
                                                  nx, ny, nz, &
                                                  dt, eps_0)
    !Input
    real, dimension(:,:,:), intent(in) :: ex_target, ey_target, ez_target
    real, dimension(:,:,:), intent(in) :: ex_source_1, ey_source_1, ez_source_1
    real, dimension(:,:,:), intent(in) :: ex_source_2, ey_source_2, ez_source_2
    real, dimension(:,:,:), intent(in) :: dx_source_1, dy_source_1, dz_source_1
    real, dimension(:,:,:), intent(in) :: dx_source_2, dy_source_2, dz_source_2
    real, dimension(:,:,:), intent(in) :: dx_source_3, dy_source_3, dz_source_3
    real, dimension(:,:,:), intent(in) :: eps_i, eps_s
    real, dimension(:,:,:), intent(in) :: tau_d, sigma
    integer, intent(in)                :: nx, ny, nz
    real, intent(in)                   :: dt, eps_0
    !type(fdtd_params_cuda), intent(in) :: params
    !type(fdtd_field_cuda), intent(in)  :: field
    !integer, value, intent(in)         :: run_num

    !Local vars
    integer                         :: ix, iy, iz
    !real, dimension(:,:,:), pointer :: ex_target
    !real, dimension(:,:,:), pointer :: ey_target
    !real, dimension(:,:,:), pointer :: ez_target

    !real, dimension(:,:,:), pointer :: ex_source_1
    !real, dimension(:,:,:), pointer :: ex_source_2
    !real, dimension(:,:,:), pointer :: dx_source_1
    !real, dimension(:,:,:), pointer :: dx_source_2
    !real, dimension(:,:,:), pointer :: dx_source_3

    !real, dimension(:,:,:), pointer :: ey_source_1
    !real, dimension(:,:,:), pointer :: ey_source_2
    !real, dimension(:,:,:), pointer :: dy_source_1
    !real, dimension(:,:,:), pointer :: dy_source_2
    !real, dimension(:,:,:), pointer :: dy_source_3

    !real, dimension(:,:,:), pointer :: ez_source_1
    !real, dimension(:,:,:), pointer :: ez_source_2
    !real, dimension(:,:,:), pointer :: dz_source_1
    !real, dimension(:,:,:), pointer :: dz_source_2
    !real, dimension(:,:,:), pointer :: dz_source_3

    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    !if(run_num .eq. 1) then
    !    ex_target => field%ex1
    !    ey_target => field%ey1
    !    ez_target => field%ez1
    !    
    !    ex_source_1 => field%ex3
    !    ex_source_2 => field%ex2
    !    dx_source_1 => field%dx1
    !    dx_source_2 => field%dx3
    !    dx_source_3 => field%dx2
    !
    !    ey_source_1 => field%ey3
    !    ey_source_2 => field%ey2
    !    dy_source_1 => field%dy1
    !    dy_source_2 => field%dy3
    !    dy_source_3 => field%dy2
    !
    !    ez_source_1 => field%ez3
    !    ez_source_2 => field%ez2
    !    dz_source_1 => field%dz1
    !    dz_source_2 => field%dz3
    !    dz_source_3 => field%dz2
    !else if(run_num .eq. 2) then
    !    ex_target => field%ex2
    !    ey_target => field%ey2
    !    ez_target => field%ez2
    !    
    !    ex_source_1 => field%ex1
    !    ex_source_2 => field%ex3
    !    dx_source_1 => field%dx2
    !    dx_source_2 => field%dx1
    !    dx_source_3 => field%dx3
    !
    !    ey_source_1 => field%ey1
    !    ey_source_2 => field%ey3
    !    dy_source_1 => field%dy2
    !    dy_source_2 => field%dy1
    !    dy_source_3 => field%dy3
    !
    !    ez_source_1 => field%ez1
    !    ez_source_2 => field%ez3
    !    dz_source_1 => field%dz2
    !    dz_source_2 => field%dz1
    !    dz_source_3 => field%dz3
    !else
    !    ex_target => field%ex3
    !    ey_target => field%ey3
    !    ez_target => field%ez3
    !    
    !    ex_source_1 => field%ex2
    !    ex_source_2 => field%ex1
    !    dx_source_1 => field%dx3
    !    dx_source_2 => field%dx2
    !    dx_source_3 => field%dx1
    !
    !    ey_source_1 => field%ey2
    !    ey_source_2 => field%ey1
    !    dy_source_1 => field%dy3
    !    dy_source_2 => field%dy2
    !    dy_source_3 => field%dy1
    !
    !    ez_source_1 => field%ez2
    !    ez_source_2 => field%ez1
    !    dz_source_1 => field%dz3
    !    dz_source_2 => field%dz2
    !    dz_source_3 => field%dz1
    !end if
    
    !Update Ex
    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz <= nz-1) then
        ex_target(ix, iy, iz) = (                                                             &
                                 1/(2 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +       &
                                 2 * dt *                                                     &
                                 (                                                            &  
                                  eps_0 * eps_s(ix, iy, iz) +                                 &
                                  sigma(ix, iy, iz) * tau_d(ix, iy, iz)                       &
                                 ) +                                                          &
                                 sigma(ix, iy, iz) * dt * dt)                                 &
                                ) *                                                           &
                                (                                                             &
                                 (                                                            &
                                  4 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +         &
                                  2 * dt *                                                    &
                                  (                                                           &
                                   eps_0 * eps_s(ix, iy, iz) +                                &
                                   sigma(ix, iy, iz) * tau_d(ix, iy, iz)                      &
                                  ) -                                                         &
                                  sigma(ix, iy, iz) * dt * dt                                 &
                                 ) *                                                          &
                                 ex_source_1(ix, iy, iz) -                                    &
                                 (2 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz)) *        &
                                 ex_source_2(ix, iy, iz) +                                    &
                                 (2 * (dt + tau_d(ix, iy, iz))) * dx_source_1(ix, iy, iz) -   &
                                 (2 * dt + 4 * tau_d(ix, iy, iz)) * dx_source_2(ix, iy, iz) + &
                                 (2*tau_d(ix, iy, iz)) * dx_source_3(ix, iy, iz)              &
                                )
    end if
    
    !Update Ey
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz <= nz-1) then
        ey_target(ix, iy, iz) = (                                                             &
                                 1/(2 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +       &
                                 2 * dt *                                                     &
                                 (                                                            &
                                  eps_0 * eps_s(ix, iy, iz) +                                 &
                                  sigma(ix, iy, iz) * tau_d(ix, iy, iz)                       &
                                 ) +                                                          &
                                 sigma(ix, iy, iz) * dt * dt)                                 &
                                ) *                                                           &
                                (                                                             &
                                 (                                                            &
                                  4 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +         &
                                  2 * dt *                                                    &  
                                  (                                                           &
                                   eps_0 * eps_s(ix, iy, iz) +                                &
                                   sigma(ix, iy, iz) * tau_d(ix, iy, iz)                      &
                                  ) -                                                         &
                                  sigma(ix, iy, iz) * dt * dt                                 &
                                 ) *                                                          &
                                 ey_source_1(ix, iy, iz) -                                    & 
                                 (2 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz)) *        &
                                 ey_source_2(ix, iy, iz) +                                    &
                                 (2 * (dt + tau_d(ix, iy, iz))) * dy_source_1(ix, iy, iz) -   &
                                 (2 * dt + 4 * tau_d(ix, iy, iz)) * dy_source_2(ix, iy, iz) + &
                                 (2 * tau_d(ix, iy, iz)) * dy_source_3(ix, iy, iz)            &
                                )
    end if
    
    !Update Ez
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz <= nz-1) then
        ez_target(ix, iy, iz) = (                                                             &
                                 1/(2 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +       &
                                 2 * dt *                                                     &
                                 (                                                            &
                                  eps_0 * eps_s(ix, iy, iz) +                                 &
                                  sigma(ix, iy, iz) * tau_d(ix, iy, iz)                       &
                                 ) +                                                          &
                                 sigma(ix, iy, iz) * dt * dt)                                 &
                                ) *                                                           &
                                (                                                             &
                                 (                                                            &
                                  4 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +         &
                                  2 * dt *                                                    &  
                                  (                                                           &
                                   eps_0 * eps_s(ix, iy, iz) +                                &
                                   sigma(ix, iy, iz) * tau_d(ix, iy, iz)                      &
                                  ) -                                                         &
                                  sigma(ix, iy, iz) * dt * dt                                 &
                                 ) *                                                          &
                                 ez_source_1(ix, iy, iz) -                                    & 
                                 (2 * eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz)) *        &
                                 ez_source_2(ix, iy, iz) +                                    &
                                 (2 * (dt + tau_d(ix, iy, iz))) * dz_source_1(ix, iy, iz) -   &
                                 (2 * dt + 4 * tau_d(ix, iy, iz)) * dz_source_2(ix, iy, iz) + &
                                 (2 * tau_d(ix, iy, iz)) * dz_source_3(ix, iy, iz)            &
                                )
    end if
end subroutine


attributes(global) subroutine update_source_cuda(dz_target, dz_source,
                                                 hx, hy,
                                                 src, jz,
                                                 dt, dx, dy, dz,
                                                 runs_count)

    !Input
    real, dimension(:,:,:), intent(in)  :: dz_target, dz_source
    real, dimension(:,:,:), intent(in)  :: hx, hy
    integer, dimension(:,:), intent(in) :: src
    real, dimension(:), intent(in)      :: jz
    real, intent(in)                    :: dt, dx, dy, dz
    integer, value, intent(in)          :: runs_count
    !type(fdtd_params_cuda), intent(in) :: params
    !type(fdtd_field_cuda), intent(in)  :: field
    !integer, value, intent(in)         :: run_num

    !Local vars 
    integer                         :: ix, iy, iz
    integer                         :: i
    integer                         :: x, y, z
    !real, dimension(:,:,:), pointer :: dz_target
    !real, dimension(:,:,:), pointer :: dz_source
    
    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    !if(run_num .eq. 1) then
    !    dz_target => field%dz1
    !    dz_source => field%dz3
    !else if(run_num .eq. 2) then
    !    dz_target => field%dz2
    !    dz_source => field%dz1
    !else
    !    dz_target => field%dz3
    !    dz_source => field%dz2
    !end if

    !Update source
    if(ix == 0 .and. iy == 0 .and. iz == 0) then
        do i=1, nsrc
            x = src(i, 1)
            y = src(i, 2)
            z = src(i, 3)
    
            dz_target(x, y, z) = dz_source(x, y, z) +                    &
                                 dt/dx * (hy(x, y, z) - hy(x-1, y, z)) - &
                                 dt/dy * (hx(x, y, z) - hx(x, y-1, z)) - &
                                 jz(((runs_count-1)*3)+1)
        end do
    end if
end subroutine


attributes(global) subroutine update_mur_boundary_cuda(ex_target, ey_target, ez_target,
                                                       ex_source, ey_source, ez_source,
                                                       rp_x_1, rp_x_end,
                                                       rp_y_1, rp_y_end,
                                                       rp_z_1, rp_z_end,
                                                       dt, dx, dy, dz,
                                                       mu_0, eps_0)

    !Input
    real, dimension(:,:,:), intent(in) :: ex_target, ey_target, ez_target
    real, dimension(:,:,:), intent(in) :: ex_source, ey_source, ez_source
    real, dimension(:,:,:), intent(in) :: rp_x_1, rp_x_end
    real, dimension(:,:,:), intent(in) :: rp_y_1, rp_y_end
    real, dimension(:,:,:), intent(in) :: rp_z_1, rp_z_end
    real, intent(in)                   :: dt, dx, dy, dz
    real, intent(in)                   :: mu_0, eps_0
    !type(fdtd_params_cuda), intent(in) :: params
    !type(fdtd_field_cuda), intent(in)  :: field
    
    !Local vars
    integer                         :: ix, iy, iz
    !real, dimension(:,:,:), pointer :: ex_target
    !real, dimension(:,:,:), pointer :: ey_target
    !real, dimension(:,:,:), pointer :: ez_target

    !real, dimension(:,:,:), pointer :: ex_source
    !real, dimension(:,:,:), pointer :: ey_source
    !real, dimension(:,:,:), pointer :: ez_source

    !Setup indexes
    ix = threadIdx%x
    iy = threadIdx%y
    iz = threadIdx%z

    !Setup based on run_num 1..3
    !if(run_num .eq. 1) then
    !    ex_target => field%ex1
    !    ey_target => field%ey1
    !    ez_target => field%ez1
    !    
    !    ex_source => field%ex3
    !    ey_source => field%ey3
    !    ez_source => field%ez3
    !else if(run_num .eq. 2) then
    !    ex_target => field%ex2
    !    ey_target => field%ey2
    !    ez_target => field%ez2
    !    
    !    ex_source => field%ex1
    !    ey_source => field%ey1
    !    ez_source => field%ez1
    !else
    !    ex_target => field%ex3
    !    ey_target => field%ey3
    !    ez_target => field%ez3
    !    
    !    ex_source => field%ex2
    !    ey_source => field%ey2
    !    ez_source => field%ez2
    !end if

    !Update Ex
    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy == 1 .and. &
       iz >= 2 .and. iz <= nz-1) then
        ex_target(ix, iy, iz) = 1/(dt + dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (dt - dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ex_target(ix, iy+1, iz) +                               &
                                 (dt + dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ex_source(ix, iy+1, iz) -                               &
                                 (dt - dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy,iz))) *    &
                                 ex_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy == ny .and. &
       iz >= 2 .and. iz <= nz-1) then
        ex_target(ix, iy, iz) = 1/(dt + dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (dt - dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                 ex_target(ix, iy-1, iz) +                                 &
                                 (dt + dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                 ex_source(ix, iy-1, iz) -                                 &
                                 (dt - dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                     &
                                )
    end if

    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz == 1) then
        ex_target(ix, iy, iz) = 1/(dt + dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ex_target(ix, iy, iz+1) +                               &
                                 (dt + dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ex_source(ix, iy, iz+1) -                               &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix >= 1 .and. ix <= nx-1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz == nz) then
        ex_target(ix, iy, iz) = 1/(dt + dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ex_target(ix, iy, iz-1) +                                 &
                                 (dt + dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ex_source(ix, iy, iz-1) -                                 &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                     &
                                )
    end if

    !Update Ey
    if(ix == 1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz <= nz-1) then
        ey_target(ix, iy, iz) = 1/(dt + dx * sqrt(mu_0 * eps_0 * rp_x_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ey_target(ix+1, iy, iz) +                               &
                                 (dt + dx * sqrt(mu_0 * eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ey_source(ix+1, iy, iz) -                               &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_1(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix == nx .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz >= 2 .and. iz <= nz-1) then
        ey_target(ix, iy, iz) = 1/(dt + dx * sqrt(mu_0 * eps_0 * rp_x_end(ix, iy, iz)))  * &
                                (                                                          &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ey_source(ix-1, iy, iz) +                                 &
                                 (dt + dx * sqrt(mu_0 * eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ey_source(ix-1, iy, iz) -                                 &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_end(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                     &
                                )
    end if

    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz == 1) then
        ey_target(ix, iy, iz) = 1/(dt + dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ey_target(ix, iy,iz+1) +                                &
                                 (dt + dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ey_source(ix, iy, iz+1) -                               &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_1(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy >= 1 .and. iy <= ny-1 .and. &
       iz == nz) then
        ey_target(ix, iy, iz) = 1/(dt + dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (dt - dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ey_target(ix, iy, iz-1) +                                 &
                                 (dt + dz * sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ey_source(ix, iy, iz-1) -                                 &
                                 (dt - dz *sqrt(mu_0 * eps_0 * rp_z_end(ix, iy, iz))) *    &
                                 ey_source(ix, iy, iz)                                     &
                                )
    end if

    !Update Ez
    if(ix == 1 .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz >= 1 .and. iz <= nz-1) then
        ez_target(ix, iy, iz) = 1/(dt + dx * sqrt(mu_0 * eps_0 * rp_x_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ez_target(ix+1, iy, iz) +                               &
                                 (dt + dx * sqrt(mu_0 * eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ez_source(ix+1, iy, iz) -                               & 
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_1(ix, iy, iz)))  *  &
                                 ez_source(ix, iy, iz)                                   &
                                )
    end if
      
    if(ix == nx .and. &
       iy >= 2 .and. iy <= ny-1 .and. &
       iz >= 1 .and. iz <= nz-1) then
        ez_target(ix, iy, iz) = 1/(dt + dx * sqrt(mu_0 * eps_0 * rp_x_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ez_target(ix-1, iy, iz) +                                 &
                                 (dt + dx * sqrt(mu_0 * eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ez_source(ix-1, iy, iz) -                                 &
                                 (dt - dx * sqrt(mu_0 * eps_0 * rp_x_end(ix, iy, iz))) *   &
                                 ez_source(ix, iy, iz)                                     &
                                )
    end if
    
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy == 1 .and. &
       iz >= 1 .and. iz <= nz-1) then 
        ez_target(ix, iy, iz) = 1/(dt + dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (dt - dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ez_target(ix, iy+1, iz) +                               & 
                                 (dt + dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ez_source(ix, iy+1, iz) -                               &
                                 (dt - dy * sqrt(mu_0 * eps_0 * rp_y_1(ix, iy, iz))) *   &
                                 ez_source(ix, iy, iz)                                   &
                                )
    end if
      
    if(ix >= 2 .and. ix <= nx-1 .and. &
       iy == ny .and. &
       iz >= 1 .and. iz <= nz-1) then 
            ez_target(ix, iy, iz) = 1/(dt + dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy, iz))) *  &
                                    (                                                          &
                                     (dt - dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                     ez_target(ix, iy-1, iz) +                                 &
                                     (dt + dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                     ez_source(ix, iy-1, iz) -                                 &
                                     (dt - dy * sqrt(mu_0 * eps_0 * rp_y_end(ix, iy, iz))) *   &
                                     ez_source(ix, iy, iz)                                     &
                                    )
    end if
end

end module
