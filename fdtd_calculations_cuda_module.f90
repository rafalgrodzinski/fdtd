module fdtd_calculations_cuda_module

use cudafor
use fdtd_data_cuda_module

implicit none

integer, parameter :: TILE_SIZE = 4
!constants
integer, constant :: g_nx, g_ny, g_nz
real, constant    :: g_dt, g_dx, g_dy, g_dz
real, constant    :: g_mu_0, g_eps_0
integer, constant :: g_nsrc

contains

attributes(global) subroutine update_h_field_cuda(hx, hy, hz,                      &
                                                  ex_source, ey_source, ez_source)

    !Input
    real, dimension(:,:,:), intent(inout) :: hx, hy, hz
    real, dimension(:,:,:), intent(in)    :: ex_source, ey_source, ez_source

    !Local vars 
    real, shared :: s_ex_source(TILE_SIZE+1, TILE_SIZE+1, TILE_SIZE+1)
    real, shared :: s_ey_source(TILE_SIZE+1, TILE_SIZE+1, TILE_SIZE+1)
    real, shared :: s_ez_source(TILE_SIZE+1, TILE_SIZE+1, TILE_SIZE+1)
    integer :: ix, iy, iz
    integer :: tix, tiy, tiz
    
    !Setup indexes
    tix = threadIdx%x
    tiy = threadIdx%y
    tiz = threadIdx%z
   
    ix = threadIdx%x + (blockIdx%x - 1) * blockDim%x
    iy = threadIdx%y + (blockIdx%y - 1) * blockDim%y
    iz = threadIdx%z + (blockIdx%z - 1) * blockDim%z
    
    !Preload data
    if(ix <= g_nx .and. iy <= g_ny .and. iz <= g_nz) then
        s_ex_source(tix, tiy, tiz) = ex_source(ix, iy, iz)
        s_ey_source(tix, tiy, tiz) = ey_source(ix, iy, iz)
        s_ez_source(tix, tiy, tiz) = ez_source(ix, iy, iz)
    end if

    if(tix == TILE_SIZE .and. ix < g_nx) then
        s_ex_source(tix+1, tiy, tiz) = ex_source(ix+1, iy, iz)
    end if
    
    if(tiy == TILE_SIZE .and. iy < g_ny) then
        s_ey_source(tix, tiy+1, tiz) = ey_source(ix, iy+1, iz)
    end if
    
    if(tiz == TILE_SIZE .and. iz < g_nz) then
        s_ez_source(tix, tiy, tiz+1) = ez_source(ix, iy, iz+1)
    end if
    
    !Wait for loads to finish
    call syncthreads()
    
    !Update Hx
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 1 .and. iz <= g_nz-1) then
        hx(ix, iy, iz) = hx(ix, iy, iz) -                                              &
                         g_dt/(g_mu_0 * g_dy) *                                              &
                         (s_ez_source(tix, tiy+1, tiz) - s_ez_source(tix, tiy, tiz)) + &
                         g_dt/(g_mu_0 * g_dz) *                                              &
                         (s_ey_source(tix, tiy, tiz+1) - s_ey_source(tix, tiy, tiz))
    end if
    
    !Update Hy
    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz >= 1 .and. iz <= g_nz-1) then
        hy(ix, iy, iz) = hy(ix, iy, iz) -                                              &
                         g_dt/(g_mu_0 * g_dz) *                                              &
                         (s_ex_source(tix, tiy, tiz+1) - s_ex_source(tix, tiy, tiz)) + &
                         g_dt/(g_mu_0 * g_dx) *                                              &
                         (s_ez_source(tix+1, tiy, tiz) - s_ez_source(tix, tiy, tiz))
    end if
    
    !Update Hz
    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        hz(ix, iy, iz) = hz(ix, iy, iz) -                                              &
                         g_dt/(g_mu_0 * g_dx) *                                              &
                         (s_ey_source(tix+1, tiy, tiz) - s_ey_source(tix, tiy, tiz)) + &
                         g_dt/(g_mu_0 * g_dy) *                                              &
                         (s_ex_source(tix, tiy+1, tiz) - s_ex_source(tix, tiy, tiz))
    end if
end subroutine


attributes(global) subroutine update_d_field_cuda(dx_target, dy_target, dz_target, &
                                                  dx_source, dy_source, dz_source, &
                                                  hx, hy, hz)

    !Input
    real, dimension(:,:,:), intent(inout) :: dx_target, dy_target, dz_target
    real, dimension(:,:,:), intent(in)    :: dx_source, dy_source, dz_source
    real, dimension(:,:,:), intent(in)    :: hx, hy, hz

    !Local vars
    integer :: ix, iy, iz
    
    !Setup indexes
    ix = threadIdx%x + (blockIdx%x - 1) * blockDim%x
    iy = threadIdx%y + (blockIdx%y - 1) * blockDim%y
    iz = threadIdx%z + (blockIdx%z - 1) * blockDim%z

    !Update g_dx
    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz<=g_nz-1) then
        dx_target(ix, iy, iz) = dx_source(ix, iy, iz) +                       &
                                g_dt/g_dy * (hz(ix, iy, iz) - hz(ix, iy-1, iz)) - &
                                g_dt/g_dz * (hy(ix, iy, iz) - hy(ix, iy, iz-1))
    end if
    
    !Update g_dy
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz<=g_nz-1) then
        dy_target(ix, iy, iz) = dy_source(ix, iy, iz) +                       &
                                g_dt/g_dz * (hx(ix, iy, iz) - hx(ix, iy, iz-1)) - &
                                g_dt/g_dx * (hz(ix, iy, iz) - hz(ix-1, iy, iz))
    end if
    
    !Update g_dz
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz >= 1 .and. iz<=g_nz-1) then
            dz_target(ix, iy, iz) = dz_source(ix, iy, iz) +                       &
                                    g_dt/g_dx * (hy(ix, iy, iz) - hy(ix-1, iy, iz)) - &
                                    g_dt/g_dy * (hx(ix, iy, iz) - hx(ix, iy-1, iz))    
    end if
end subroutine


attributes(global) subroutine update_e_field_cuda(ex_target, ey_target, ez_target,       &
                                                  ex_source_1, ey_source_1, ez_source_1, &
                                                  ex_source_2, ey_source_2, ez_source_2, &
                                                  dx_source_1, dy_source_1, dz_source_1, &
                                                  dx_source_2, dy_source_2, dz_source_2, &
                                                  dx_source_3, dy_source_3, dz_source_3, &
                                                  eps_i, eps_s,                          &
                                                  tau_d, sigma)

    !Input
    real, dimension(:,:,:), intent(inout) :: ex_target, ey_target, ez_target
    real, dimension(:,:,:), intent(in)    :: ex_source_1, ey_source_1, ez_source_1
    real, dimension(:,:,:), intent(in)    :: ex_source_2, ey_source_2, ez_source_2
    real, dimension(:,:,:), intent(in)    :: dx_source_1, dy_source_1, dz_source_1
    real, dimension(:,:,:), intent(in)    :: dx_source_2, dy_source_2, dz_source_2
    real, dimension(:,:,:), intent(in)    :: dx_source_3, dy_source_3, dz_source_3
    real, dimension(:,:,:), intent(in)    :: eps_i, eps_s
    real, dimension(:,:,:), intent(in)    :: tau_d, sigma
    !Local vars
    integer :: ix, iy, iz

    !Setup indexes
    ix = threadIdx%x + (blockIdx%x - 1) * blockDim%x
    iy = threadIdx%y + (blockIdx%y - 1) * blockDim%y
    iz = threadIdx%z + (blockIdx%z - 1) * blockDim%z

    !Update Ex
    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ex_target(ix, iy, iz) = (                                                             &
                                 1/(2 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +       &
                                 2 * g_dt *                                                     &
                                 (                                                            &  
                                  g_eps_0 * eps_s(ix, iy, iz) +                                 &
                                  sigma(ix, iy, iz) * tau_d(ix, iy, iz)                       &
                                 ) +                                                          &
                                 sigma(ix, iy, iz) * g_dt * g_dt)                                 &
                                ) *                                                           &
                                (                                                             &
                                 (                                                            &
                                  4 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +         &
                                  2 * g_dt *                                                    &
                                  (                                                           &
                                   g_eps_0 * eps_s(ix, iy, iz) +                                &
                                   sigma(ix, iy, iz) * tau_d(ix, iy, iz)                      &
                                  ) -                                                         &
                                  sigma(ix, iy, iz) * g_dt * g_dt                                 &
                                 ) *                                                          &
                                 ex_source_1(ix, iy, iz) -                                    &
                                 (2 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz)) *        &
                                 ex_source_2(ix, iy, iz) +                                    &
                                 (2 * (g_dt + tau_d(ix, iy, iz))) * dx_source_1(ix, iy, iz) -   &
                                 (2 * g_dt + 4 * tau_d(ix, iy, iz)) * dx_source_2(ix, iy, iz) + &
                                 (2*tau_d(ix, iy, iz)) * dx_source_3(ix, iy, iz)              &
                                )
    end if
    
    !Update Ey
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ey_target(ix, iy, iz) = (                                                             &
                                 1/(2 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +       &
                                 2 * g_dt *                                                     &
                                 (                                                            &
                                  g_eps_0 * eps_s(ix, iy, iz) +                                 &
                                  sigma(ix, iy, iz) * tau_d(ix, iy, iz)                       &
                                 ) +                                                          &
                                 sigma(ix, iy, iz) * g_dt * g_dt)                                 &
                                ) *                                                           &
                                (                                                             &
                                 (                                                            &
                                  4 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +         &
                                  2 * g_dt *                                                    &  
                                  (                                                           &
                                   g_eps_0 * eps_s(ix, iy, iz) +                                &
                                   sigma(ix, iy, iz) * tau_d(ix, iy, iz)                      &
                                  ) -                                                         &
                                  sigma(ix, iy, iz) * g_dt * g_dt                                 &
                                 ) *                                                          &
                                 ey_source_1(ix, iy, iz) -                                    & 
                                 (2 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz)) *        &
                                 ey_source_2(ix, iy, iz) +                                    &
                                 (2 * (g_dt + tau_d(ix, iy, iz))) * dy_source_1(ix, iy, iz) -   &
                                 (2 * g_dt + 4 * tau_d(ix, iy, iz)) * dy_source_2(ix, iy, iz) + &
                                 (2 * tau_d(ix, iy, iz)) * dy_source_3(ix, iy, iz)            &
                                )
    end if
    
    !Update Ez
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ez_target(ix, iy, iz) = (                                                             &
                                 1/(2 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +       &
                                 2 * g_dt *                                                     &
                                 (                                                            &
                                  g_eps_0 * eps_s(ix, iy, iz) +                                 &
                                  sigma(ix, iy, iz) * tau_d(ix, iy, iz)                       &
                                 ) +                                                          &
                                 sigma(ix, iy, iz) * g_dt * g_dt)                                 &
                                ) *                                                           &
                                (                                                             &
                                 (                                                            &
                                  4 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz) +         &
                                  2 * g_dt *                                                    &  
                                  (                                                           &
                                   g_eps_0 * eps_s(ix, iy, iz) +                                &
                                   sigma(ix, iy, iz) * tau_d(ix, iy, iz)                      &
                                  ) -                                                         &
                                  sigma(ix, iy, iz) * g_dt * g_dt                                 &
                                 ) *                                                          &
                                 ez_source_1(ix, iy, iz) -                                    & 
                                 (2 * g_eps_0 * eps_i(ix, iy, iz) * tau_d(ix, iy, iz)) *        &
                                 ez_source_2(ix, iy, iz) +                                    &
                                 (2 * (g_dt + tau_d(ix, iy, iz))) * dz_source_1(ix, iy, iz) -   &
                                 (2 * g_dt + 4 * tau_d(ix, iy, iz)) * dz_source_2(ix, iy, iz) + &
                                 (2 * tau_d(ix, iy, iz)) * dz_source_3(ix, iy, iz)            &
                                )
    end if
end subroutine


attributes(global) subroutine update_source_cuda(dz_target, dz_source, &
                                                 hx, hy,               &
                                                 src, jz,              &
                                                 runs_count)

    !Input
    real, dimension(:,:,:), intent(inout) :: dz_target
    real, dimension(:,:,:), intent(in)    :: dz_source
    real, dimension(:,:,:), intent(in)    :: hx, hy
    integer, dimension(:,:), intent(in)   :: src
    real, dimension(:), intent(in)        :: jz
    integer, value, intent(in)            :: runs_count

    !Local vars 
    integer :: ix, iy, iz
    integer :: i
    integer :: x, y, z
    
    !Setup indexes
    ix = threadIdx%x + (blockIdx%x - 1) * blockDim%x
    iy = threadIdx%y + (blockIdx%y - 1) * blockDim%y
    iz = threadIdx%z + (blockIdx%z - 1) * blockDim%z

    !Update source
    if(ix == 1 .and. iy == 1 .and. iz == 1) then
        do i=1, g_nsrc_0
            x = src(i, 1)
            y = src(i, 2)
            z = src(i, 3)
    
            dz_target(x, y, z) = dz_source(x, y, z) +                    &
                                 g_dt/g_dx * (hy(x, y, z) - hy(x-1, y, z)) - &
                                 g_dt/g_dy * (hx(x, y, z) - hx(x, y-1, z)) - &
                                 jz(((runs_count-1)*3)+1)
        end do
    end if
end subroutine


attributes(global) subroutine update_mur_boundary_cuda(ex_target, ey_target, ez_target, &
                                                       ex_source, ey_source, ez_source, &
                                                       rp_x_1, rp_x_end,                &
                                                       rp_y_1, rp_y_end,                &
                                                       rp_z_1, rp_z_end)

    !Input
    real, dimension(:,:,:), intent(inout) :: ex_target, ey_target, ez_target
    real, dimension(:,:,:), intent(in)    :: ex_source, ey_source, ez_source
    real, dimension(:,:,:), intent(in)    :: rp_x_1, rp_x_end
    real, dimension(:,:,:), intent(in)    :: rp_y_1, rp_y_end
    real, dimension(:,:,:), intent(in)    :: rp_z_1, rp_z_end
    
    !Local vars
    integer :: ix, iy, iz

    !Setup indexes
    ix = threadIdx%x + (blockIdx%x - 1) * blockDim%x
    iy = threadIdx%y + (blockIdx%y - 1) * blockDim%y
    iz = threadIdx%z + (blockIdx%z - 1) * blockDim%z

    !Update Ex
    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy == 1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ex_target(ix, iy, iz) = 1/(g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ex_target(ix, iy+1, iz) +                               &
                                 (g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ex_source(ix, iy+1, iz) -                               &
                                 (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy,iz))) *    &
                                 ex_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy == g_ny .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ex_target(ix, iy, iz) = 1/(g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                 ex_target(ix, iy-1, iz) +                                 &
                                 (g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                 ex_source(ix, iy-1, iz) -                                 &
                                 (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                     &
                                )
    end if

    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz == 1) then
        ex_target(ix, iy, iz) = 1/(g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ex_target(ix, iy, iz+1) +                               &
                                 (g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ex_source(ix, iy, iz+1) -                               &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix >= 1 .and. ix <= g_nx-1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz == g_nz) then
        ex_target(ix, iy, iz) = 1/(g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ex_target(ix, iy, iz-1) +                                 &
                                 (g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ex_source(ix, iy, iz-1) -                                 &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz))) *   &
                                 ex_source(ix, iy, iz)                                     &
                                )
    end if

    !Update Ey
    if(ix == 1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ey_target(ix, iy, iz) = 1/(g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ey_target(ix+1, iy, iz) +                               &
                                 (g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ey_source(ix+1, iy, iz) -                               &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix == g_nx .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz >= 2 .and. iz <= g_nz-1) then
        ey_target(ix, iy, iz) = 1/(g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix, iy, iz)))  * &
                                (                                                          &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ey_source(ix-1, iy, iz) +                                 &
                                 (g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ey_source(ix-1, iy, iz) -                                 &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                     &
                                )
    end if

    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz == 1) then
        ey_target(ix, iy, iz) = 1/(g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ey_target(ix, iy,iz+1) +                                &
                                 (g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz+1))) * &
                                 ey_source(ix, iy, iz+1) -                               &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_1(ix, iy, iz))) *   &
                                 ey_source(ix, iy, iz)                                   &
                                )
    end if

    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy >= 1 .and. iy <= g_ny-1 .and. &
       iz == g_nz) then
        ey_target(ix, iy, iz) = 1/(g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (g_dt - g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ey_target(ix, iy, iz-1) +                                 &
                                 (g_dt + g_dz * sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz-1))) * &
                                 ey_source(ix, iy, iz-1) -                                 &
                                 (g_dt - g_dz *sqrt(g_mu_0 * g_eps_0 * rp_z_end(ix, iy, iz))) *    &
                                 ey_source(ix, iy, iz)                                     &
                                )
    end if

    !Update Ez
    if(ix == 1 .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz >= 1 .and. iz <= g_nz-1) then
        ez_target(ix, iy, iz) = 1/(g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ez_target(ix+1, iy, iz) +                               &
                                 (g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix+1, iy, iz))) * &
                                 ez_source(ix+1, iy, iz) -                               & 
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_1(ix, iy, iz)))  *  &
                                 ez_source(ix, iy, iz)                                   &
                                )
    end if
      
    if(ix == g_nx .and. &
       iy >= 2 .and. iy <= g_ny-1 .and. &
       iz >= 1 .and. iz <= g_nz-1) then
        ez_target(ix, iy, iz) = 1/(g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix, iy, iz))) *  &
                                (                                                          &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ez_target(ix-1, iy, iz) +                                 &
                                 (g_dt + g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix-1, iy, iz))) * &
                                 ez_source(ix-1, iy, iz) -                                 &
                                 (g_dt - g_dx * sqrt(g_mu_0 * g_eps_0 * rp_x_end(ix, iy, iz))) *   &
                                 ez_source(ix, iy, iz)                                     &
                                )
    end if
    
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy == 1 .and. &
       iz >= 1 .and. iz <= g_nz-1) then 
        ez_target(ix, iy, iz) = 1/(g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy, iz))) *  &
                                (                                                        &
                                 (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ez_target(ix, iy+1, iz) +                               & 
                                 (g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy+1, iz))) * &
                                 ez_source(ix, iy+1, iz) -                               &
                                 (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_1(ix, iy, iz))) *   &
                                 ez_source(ix, iy, iz)                                   &
                                )
    end if
      
    if(ix >= 2 .and. ix <= g_nx-1 .and. &
       iy == g_ny .and. &
       iz >= 1 .and. iz <= g_nz-1) then 
            ez_target(ix, iy, iz) = 1/(g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy, iz))) *  &
                                    (                                                          &
                                     (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                     ez_target(ix, iy-1, iz) +                                 &
                                     (g_dt + g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy-1, iz))) * &
                                     ez_source(ix, iy-1, iz) -                                 &
                                     (g_dt - g_dy * sqrt(g_mu_0 * g_eps_0 * rp_y_end(ix, iy, iz))) *   &
                                     ez_source(ix, iy, iz)                                     &
                                    )
    end if
end subroutine

end module
