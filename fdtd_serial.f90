include "fdtd_data_module.f90"


program fdtd_serial
use fdtd_data_module

implicit none

type(fdtd_state), pointer :: state
type(fdtd_field), pointer :: field

call init_fdtd_state(state, "abc")

end program