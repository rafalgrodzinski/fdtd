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
    
    !setup based on run_num 0..2
    if(run_num .eq. 0) then
        ex_field => field%ex3
        ey_field => field%ey3
        ez_field => field%ez3
    else if(run_num .eq. 1) then
        ex_field => field%ex1
        ey_field => field%ey1
        ez_field => field%ez1
    else
        ex_field => field%ex2
        ey_field => field%ey2
        ez_field => field%ez2
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
    real, dimension(:,:,:), pointer :: dx_target
    real, dimension(:,:,:), pointer :: dy_target
    real, dimension(:,:,:), pointer :: dz_target
    
    real, dimension(:,:,:), pointer :: dx_source
    real, dimension(:,:,:), pointer :: dy_source
    real, dimension(:,:,:), pointer :: dz_source
    
    !setup based on run_num 0..2
    if(run_num .eq. 0) then
        dx_target => field%dx1
        dy_target => field%dx1
        dz_target => field%dx1
        
        dx_source => field%dx3
        dy_source => field%dx3
        dz_source => field%dx3
    else if(run_num .eq. 1) then
        dx_target => field%dx2
        dy_target => field%dx2
        dz_target => field%dx2
        
        dx_source => field%dx1
        dy_source => field%dx1        
        dz_source => field%dx1
    else
        dx_target => field%dx3
        dy_target => field%dx3
        dz_target => field%dx3

        dx_source => field%dx2        
        dy_source => field%dx2        
        dz_source => field%dx2
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


subroutine update_e_field(field, state, run_num)
    !input
    type(fdtd_field), pointer, intent(inout) :: field
    type(fdtd_state), pointer, intent(in)    :: state
    integer, intent(in)                      :: run_num
    !local vars
    integer :: ix, iy, iz
    real, dimension(:,:,:), pointer :: ex_target
    real, dimension(:,:,:), pointer :: ey_target
    real, dimension(:,:,:), pointer :: ez_target
    
    real, dimension(:,:,:), pointer :: ex_source_1
    real, dimension(:,:,:), pointer :: ex_source_2
    real, dimension(:,:,:), pointer :: dx_source_1
    real, dimension(:,:,:), pointer :: dx_source_2
    real, dimension(:,:,:), pointer :: dx_source_3
    
    real, dimension(:,:,:), pointer :: ey_source_1
    real, dimension(:,:,:), pointer :: ey_source_2
    real, dimension(:,:,:), pointer :: dy_source_1
    real, dimension(:,:,:), pointer :: dy_source_2
    real, dimension(:,:,:), pointer :: dy_source_3
    
    real, dimension(:,:,:), pointer :: ez_source_1
    real, dimension(:,:,:), pointer :: ez_source_2
    real, dimension(:,:,:), pointer :: dz_source_1
    real, dimension(:,:,:), pointer :: dz_source_2
    real, dimension(:,:,:), pointer :: dz_source_3
    
    !setup based on run_num 0..2
    if(run_num .eq. 0) then
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
    else if(run_num .eq. 1) then
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
    do iz=2, state%nz-1
        do iy=2, state%ny-1
	        do ix=1, state%nx-1
                ex_target(ix, iy, iz) = (                                                                         &
                                         1/(2 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) + &
                                         2 * state%dt *                                                           &
                                         (                                                                        &  
                                          state%eps_0 * field%eps_s(ix, iy, iz) +                                 &
     		                              field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                       &
     		                             ) +                                                                      &
     	                                 field%sigma(ix, iy, iz) * state%dt * state%dt)                           &
     	                                ) *                                                                       &
     	                                (                                                                         &
     	                                 (                                                                        &
     	                                  4 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) +   &
     	                                  2 * state%dt *                                                          &
     	                                  (                                                                       &
     	                                   state%eps_0 * field%eps_s(ix, iy, iz) +                                &
     	                                   field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                      &
     	                                  ) -                                                                     &
      	                                  field%sigma(ix, iy, iz) * state%dt * state%dt                           &
      	                                 ) *                                                                      &
                                         ex_source_1(ix, iy, iz) -                                                &
     		                             (2 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz)) *  & 
                                         ex_source_2(ix, iy, iz) +                                                &
     		                             (2 * (state%dt + field%tau_d(ix, iy, iz))) * dx_source_1(ix, iy, iz) -   &
     	                                 (2 * state%dt + 4 * field%tau_d(ix, iy, iz)) * dx_source_2(ix, iy, iz) + &
     		                             (2*field%tau_d(ix, iy, iz)) * dx_source_3(ix, iy, iz)                    &
     		                            )
            end do
	    end do
    end do
    
    !Update Ey
    do iz=2, state%nz-1
        do iy=1, state%ny-1
	        do ix=2, state%nx-1
                ey_target(ix, iy, iz) = (                                                                         &
                                         1/(2 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) + &
                                         2 * state%dt *                                                           &
                                         (                                                                        &
                                          state%eps_0 * field%eps_s(ix, iy, iz) +                                 &
     		                              field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                       &
     		                             ) +                                                                      &
     	                                 field%sigma(ix, iy, iz) * state%dt * state%dt)                           &
     	                                ) *                                                                       &
                                        (                                                                         &
                                         (                                                                        &
                                          4 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) +   &
                    		              2 * state%dt *                                                          &  
                    		              (                                                                       &
                    		               state%eps_0 * field%eps_s(ix, iy, iz) +                                &
                    		               field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                      &
                    		              ) -                                                                     &
                    		              field%sigma(ix, iy, iz) * state%dt * state%dt                           &
                    		             ) *                                                                      &
                    		             ey_source_1(ix, iy, iz) -                                                & 
                                		 (2 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz)) *  &
                                         ey_source_2(ix, iy, iz) +                                                &
                                		 (2 * (state%dt + field%tau_d(ix, iy, iz))) * dy_source_1(ix, iy, iz) -   &
     		                             (2 * state%dt + 4 * field%tau_d(ix, iy, iz)) * dy_source_2(ix, iy, iz) + &
                                		 (2 * field%tau_d(ix, iy, iz)) * dy_source_3(ix, iy, iz)                  &
                                		)
            end do
	    end do
    end do
    
    !Update Ez
    do iz=2, state%nz-1
        do iy=1, state%ny-1
	        do ix=2, state%nx-1
                ez_target(ix, iy, iz) = (                                                                         &
                                         1/(2 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) + &
                                         2 * state%dt *                                                           &
                                         (                                                                        &
                                          state%eps_0 * field%eps_s(ix, iy, iz) +                                 &
     		                              field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                       &
     		                             ) +                                                                      &
     	                                 field%sigma(ix, iy, iz) * state%dt * state%dt)                           &
     	                                ) *                                                                       &
                                        (                                                                         &
                                         (                                                                        &
                                          4 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz) +   &
                    		              2 * state%dt *                                                          &  
                    		              (                                                                       &
                    		               state%eps_0 * field%eps_s(ix, iy, iz) +                                &
                    		               field%sigma(ix, iy, iz) * field%tau_d(ix, iy, iz)                      &
                    		              ) -                                                                     &
                    		              field%sigma(ix, iy, iz) * state%dt * state%dt                           &
                    		             ) *                                                                      &
                    		             ez_source_1(ix, iy, iz) -                                                & 
                                		 (2 * state%eps_0 * field%eps_i(ix, iy, iz) * field%tau_d(ix, iy, iz)) *  &
                                         ez_source_2(ix, iy, iz) +                                                &
                                		 (2 * (state%dt + field%tau_d(ix, iy, iz))) * dz_source_1(ix, iy, iz) -   &
     		                             (2 * state%dt + 4 * field%tau_d(ix, iy, iz)) * dz_source_2(ix, iy, iz) + &
                                		 (2 * field%tau_d(ix, iy, iz)) * dz_source_3(ix, iy, iz)                  &
                                		)
            end do
	    end do
    end do
end

end