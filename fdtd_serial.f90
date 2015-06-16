include "utils_module.f90"
include "fdtd_data_module.f90"
include "fdtd_calculations_module.f90"


program fdtd_serial
use fdtd_data_module
use fdtd_calculations_module

implicit none

type(fdtd_state), pointer :: state
type(fdtd_field), pointer :: field

call init_fdtd_state(state, "data/input_params")
call init_fdtd_field(field, state)

end program