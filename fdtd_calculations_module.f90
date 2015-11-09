module fdtd_calculations_module

use fdtd_data_module


implicit none

contains

subroutine update_h_field(params, field,                   &
						  ex_source, ey_source, ez_source, &
						  run_num)

    !Input
    type(fdtd_params), intent(inout)   :: params
    type(fdtd_field), intent(inout)    :: field

	real, dimension(:,:,:), intent(in) :: ex_source
    real, dimension(:,:,:), intent(in) :: ey_source
    real, dimension(:,:,:), intent(in) :: ez_source

    integer, intent(in)                :: run_num

    !Local vars
    integer :: ix, iy, iz
    
    !Update Hx
    do iz=1, params%nz-1
        do iy=1, params%ny-1
	        do ix=2, params%nx-1
                field%hx(ix, iy, iz) = field%hx(ix, iy, iz) -                              &
                                       params%dt/(params%mu_0 * params%dy) *               &
                                       (ez_source(ix, iy+1, iz) - ez_source(ix, iy, iz)) + &
                                       params%dt/(params%mu_0 * params%dz) *               &
                                       (ey_source(ix, iy, iz+1) - ey_source(ix, iy, iz))
            end do
	    end do
    end do

    !Update Hy
    do iz=1, params%nz-1
        do iy=2, params%ny-1
	        do ix=1, params%nx-1
                field%hy(ix, iy, iz) = field%hy(ix, iy, iz) -                              &
                                       params%dt/(params%mu_0 * params%dz) *               &
                                       (ex_source(ix, iy, iz+1) - ex_source(ix, iy, iz)) + &
                                       params%dt/(params%mu_0 * params%dx) *               &
                                       (ez_source(ix+1, iy, iz) - ez_source(ix, iy, iz))
            end do
	    end do
    end do

    !Update Hz
    do iz=2, params%nz-1
        do iy=1, params%ny-1
	        do ix=1, params%nx-1
                field%hz(ix, iy, iz) = field%hz(ix, iy, iz) -                              &
                                       params%dt/(params%mu_0 * params%dx) *               &
                                       (ey_source(ix+1, iy, iz) - ey_source(ix, iy, iz)) + &
                                       params%dt/(params%mu_0 * params%dy) *               &
                                       (ex_source(ix, iy+1, iz) - ex_source(ix, iy, iz))
            end do
	    end do
    end do
end subroutine


subroutine update_d_field(params, field,                   &
 						  dx_target, dy_target, dz_target, &
 						  dx_source, dy_source, dz_source, &
						  run_num)

    !Input
    type(fdtd_params), intent(inout) 	  :: params
    type(fdtd_field), intent(inout)  	  :: field

	real, dimension(:,:,:), intent(inout) :: dx_target
    real, dimension(:,:,:), intent(inout) :: dy_target
    real, dimension(:,:,:), intent(inout) :: dz_target

    real, dimension(:,:,:), intent(in)    :: dx_source
    real, dimension(:,:,:), intent(in)    :: dy_source
    real, dimension(:,:,:), intent(in)    :: dz_source

    integer, intent(in)              	  :: run_num

    !Local vars
    integer :: ix, iy, iz

    !Update Dx
    do iz=2, params%nz-1
        do iy=2, params%ny-1
	        do ix=1, params%nx-1
                dx_target(ix, iy, iz) = dx_source(ix, iy, iz) +                                                 &
                                        params%dt/params%dy * (field%hz(ix, iy, iz) - field%hz(ix, iy-1, iz)) - &
                                        params%dt/params%dz * (field%hy(ix, iy, iz) - field%hy(ix, iy, iz-1))
            end do
	    end do
    end do

    !Update Dy
    do iz=2, params%nz-1
        do iy=1, params%ny-1
	        do ix=2, params%nx-1
                dy_target(ix, iy, iz) = dy_source(ix, iy, iz) +                                                 &
                                        params%dt/params%dz * (field%hx(ix, iy, iz) - field%hx(ix, iy, iz-1)) - &
                                        params%dt/params%dx * (field%hz(ix, iy, iz) - field%hz(ix-1, iy, iz))
            end do
	    end do
    end do

    !Update Dz
    do iz=1, params%nz-1
        do iy=2, params%ny-1
            do ix=2, params%nx-1
                dz_target(ix, iy, iz) = dz_source(ix, iy, iz) +                                                 &
                                        params%dt/params%dx * (field%hy(ix, iy, iz) - field%hy(ix-1, iy, iz)) - &
                                        params%dt/params%dy * (field%hx(ix, iy, iz) - field%hx(ix, iy-1, iz))                                
            end do
        end do
    end do
end subroutine


subroutine update_e_field(params, field,                         &
                          ex_target, ey_target, ez_target,       &
                          ex_source_1, ey_source_1, ez_source_1, &
                          ex_source_2, ey_source_2, ez_source_2, &
                          dx_source_1, dy_source_1, dz_source_1, &
                          dx_source_2, dy_source_2, dz_source_2, &
                          dx_source_3, dy_source_3, dz_source_3, &
                          run_num)

    !Input
    type(fdtd_params), intent(inout)      :: params
    type(fdtd_field), intent(inout)       :: field

	real, dimension(:,:,:), intent(inout) :: ex_target
    real, dimension(:,:,:), intent(inout) :: ey_target
    real, dimension(:,:,:), intent(inout) :: ez_target

	real, dimension(:,:,:), intent(in)    :: ex_source_1
	real, dimension(:,:,:), intent(in)    :: ey_source_1
	real, dimension(:,:,:), intent(in)    :: ez_source_1
	
	real, dimension(:,:,:), intent(in)    :: ex_source_2
	real, dimension(:,:,:), intent(in)    :: ey_source_2
	real, dimension(:,:,:), intent(in)    :: ez_source_2

    real, dimension(:,:,:), intent(in)    :: dx_source_1
    real, dimension(:,:,:), intent(in)    :: dy_source_1
    real, dimension(:,:,:), intent(in)    :: dz_source_1
    
    real, dimension(:,:,:), intent(in)    :: dx_source_2
    real, dimension(:,:,:), intent(in)    :: dy_source_2
    real, dimension(:,:,:), intent(in)    :: dz_source_2
    
    real, dimension(:,:,:), intent(in)    :: dx_source_3
    real, dimension(:,:,:), intent(in)    :: dy_source_3
    real, dimension(:,:,:), intent(in)    :: dz_source_3
    
	integer, intent(in) 			      :: run_num

    !Local vars
    integer :: ix, iy, iz

    !Update Ex
    do iz=2, params%nz-1
        do iy=2, params%ny-1
	        do ix=1, params%nx-1
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
            end do
	    end do
    end do

    !Update Ey
    do iz=2, params%nz-1
        do iy=1, params%ny-1
	        do ix=2, params%nx-1
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
            end do
	    end do
    end do

    !Update Ez
    do iz=2, params%nz-1
        do iy=1, params%ny-1
	        do ix=2, params%nx-1
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
            end do
	    end do
    end do
end subroutine


subroutine update_source(params, field,        &
						 dz_target, dz_source, &
                         run_num, runs_count)

    !Input
    type(fdtd_params), intent(inout)      :: params
    type(fdtd_field), intent(inout)       :: field
    
	real, dimension(:,:,:), intent(inout) :: dz_target
    real, dimension(:,:,:), intent(in)    :: dz_source
    
    integer, intent(in)                   :: run_num
    integer, intent(in)                   :: runs_count

    !Local vars
    integer :: i
    integer :: x, y, z

    !Update source
    do i=1, params%nsrc
        x = params%src(i, 1)
        y = params%src(i, 2)
        z = params%src(i, 3)
    
        dz_target(x, y, z) = dz_source(x, y, z) +                     &
            params%dt/params%dx * (field%hy(x, y, z) - field%hy(x-1, y, z)) - &
            params%dt/params%dy * (field%hx(x, y, z) - field%hx(x, y-1, z)) - &
            params%jz(((runs_count-1)*3)+1)
    end do
end subroutine


subroutine update_mur_boundary(params, field,                   &
                               ex_target, ey_target, ez_target, &
                               ex_source, ey_source, ez_source, &
                               run_num)

    !Input
    type(fdtd_params), intent(inout)      :: params
    type(fdtd_field), intent(inout)       :: field

	real, dimension(:,:,:), intent(inout) :: ex_target
    real, dimension(:,:,:), intent(inout) :: ey_target
    real, dimension(:,:,:), intent(inout) :: ez_target

    real, dimension(:,:,:), intent(in)    :: ex_source
    real, dimension(:,:,:), intent(in)    :: ey_source
    real, dimension(:,:,:), intent(in)    :: ez_source

    integer, intent(in)                   :: run_num

    !Local vars
    integer :: ix, iy, iz

    !Update Ex
    iy=1
    do iz=2, params%nz-1
        do ix=1, params%nx-1
            ex_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy, iz))) *  &
                                    (                                                                                      &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                     ex_target(ix, iy+1, iz) +                                                             &
                                     (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                     ex_source(ix, iy+1, iz) -                                                             &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy,iz))) *    &
                                     ex_source(ix, iy, iz)                                                                 &
                                    )
        end do
    end do
    
    iy=params%ny
    do iz=2, params%nz-1
        do ix=1, params%nx-1
            ex_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *  &
                                    (                                                                                        &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                     ex_target(ix, iy-1, iz) +                                                               &
                                     (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                     ex_source(ix, iy-1, iz) -                                                               &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *   &
                                     ex_source(ix, iy, iz)                                                                   &
                                    )
        end do
    end do
    
    iz=1
    do iy=2, params%ny-1
        do ix=1, params%nx-1
            ex_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *  &
                                    (                                                                                      &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                     ex_target(ix, iy, iz+1) +                                                             &
                                     (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                     ex_source(ix, iy, iz+1) -                                                             &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *   &
                                     ex_source(ix, iy, iz)                                                                 &
                                    )
        end do
    end do
    
    iz=params%nz
    do iy=2, params%ny-1
        do ix=1, params%nx-1
            ex_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *  &
                                    (                                                                                        &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                     ex_target(ix, iy, iz-1) +                                                               &
                                     (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                     ex_source(ix, iy, iz-1) -                                                               &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *   & 
                                     ex_source(ix, iy, iz)                                                                   &
                                    )
        end do
    end do
    
    !Update Ey
    ix=1
    do iz=2, params%nz-1
        do iy=1, params%ny-1
            ey_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz))) *  &
                                    (                                                                                      &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * &
                                     ey_target(ix+1, iy, iz) +                                                             &
                                     (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * &
                                     ey_source(ix+1, iy, iz) -                                                             &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz))) *   &
                                     ey_source(ix, iy, iz)                                                                 &
                                    )
        end do
    end do
      
    ix=params%nx
    do iz=2, params%nz-1
        do iy=1, params%ny-1
            ey_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz)))  * &
                                    (                                                                                        &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                     ey_source(ix-1, iy, iz) +                                                               &
                                     (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                     ey_source(ix-1, iy, iz) -                                                               &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz))) *   &
                                     ey_source(ix, iy, iz)                                                                   &
                                    )
        end do
    end do
      
    iz=1
    do iy=1, params%ny-1
        do ix=2, params%nx-1
            ey_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *  &
                                    (                                                                                      &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                     ey_target(ix, iy,iz+1) +                                                              &
                                     (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz+1))) * &
                                     ey_source(ix, iy, iz+1) -                                                             &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_1(ix, iy, iz))) *   &
                                     ey_source(ix, iy, iz)                                                                 &
                                    )
        end do
    end do
      
    iz=params%nz
    do iy=1, params%ny-1
        do ix=2, params%nx-1
            ey_target(ix, iy, iz) = 1/(params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *  &
                                    (                                                                                        &
                                     (params%dt - params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                     ey_target(ix, iy, iz-1) +                                                               &
                                     (params%dt + params%dz * sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz-1))) * &
                                     ey_source(ix, iy, iz-1) -                                                               &
                                     (params%dt - params%dz *sqrt(params%mu_0 * params%eps_0 * field%rp_z_end(ix, iy, iz))) *    &     
                                     ey_source(ix, iy, iz)                                                                   &
                                    )
        end do
    end do
    
    !Update Ez
    ix=1
    do iz=1, params%nz-1
        do iy=2, params%ny-1
            ez_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz))) *  &
                                    (                                                                                      &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * & 
                                     ez_target(ix+1, iy, iz) +                                                             &
                                     (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix+1, iy, iz))) * &
                                     ez_source(ix+1, iy, iz) -                                                             & 
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_1(ix, iy, iz)))  *  &
                                     ez_source(ix, iy, iz)                                                                 &
                                    )

        end do
    end do
      
    ix=params%nx
    do iz=1, params%nz-1
        do iy=2, params%ny-1
            ez_target(ix, iy, iz) = 1/(params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz))) *  &
                                    (                                                                                        &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                     ez_target(ix-1, iy, iz) +                                                               &
                                     (params%dt + params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix-1, iy, iz))) * &
                                     ez_source(ix-1, iy, iz) -                                                               &
                                     (params%dt - params%dx * sqrt(params%mu_0 * params%eps_0 * field%rp_x_end(ix, iy, iz))) *   &
                                     ez_source(ix, iy, iz)                                                                   &
                                    )
         end do
      end do
      
    iy=1
    do iz=1, params%nz-1
        do ix=2, params%nx-1
            ez_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy, iz))) *  &
                                    (                                                                                      &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                     ez_target(ix, iy+1, iz) +                                                             & 
                                     (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy+1, iz))) * &
                                     ez_source(ix, iy+1, iz) -                                                             &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_1(ix, iy, iz))) *   &
                                     ez_source(ix, iy, iz)                                                                 &
                                    )
        end do
    end do
      
    iy=params%ny
    do iz=1, params%nz-1
        do ix=2, params%nx-1
            ez_target(ix, iy, iz) = 1/(params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *  &
                                    (                                                                                        &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                     ez_target(ix, iy-1, iz) +                                                               &
                                     (params%dt + params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy-1, iz))) * &
                                     ez_source(ix, iy-1, iz) -                                                               &
                                     (params%dt - params%dy * sqrt(params%mu_0 * params%eps_0 * field%rp_y_end(ix, iy, iz))) *   &
                                     ez_source(ix, iy, iz)                                                                   &
                                    )
        end do
    end do
end subroutine

end module
