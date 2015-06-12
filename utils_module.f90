module utils_module
implicit none

contains

subroutine check_error(error_code, message)
    integer, intent(in)          :: error_code
    character(len=*), intent(in) :: message

    if(error_code .ne. 0) then
        print '(A, I3, A, A)', "Error, Code ", error_code, ": ", message
        stop
    endif

end

end