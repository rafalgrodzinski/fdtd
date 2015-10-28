module utils_module

implicit none

contains

subroutine check_error(error_code, message)
    !inut
    integer, intent(in)          :: error_code
    character(len=*), intent(in) :: message

    if(error_code .ne. 0) then
        print '(A, I3, A, A)', "Error, Code ", error_code, ": ", message
        stop
    endif
end subroutine


character(len=128) function str(i)
    !input
    integer, intent(in) :: i
    
    write (str, *) i
    str = adjustl(str)
end function


character(len=128) function generate_file_name(prefix, postfix, number)
    !input
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: postfix
    integer, intent(in)          :: number
    
    write(generate_file_name, fmt='(I5)'), number
    generate_file_name = "0000" // adjustl(generate_file_name)
    generate_file_name = generate_file_name(len(trim(generate_file_name))-4 : len(trim(generate_file_name)))
    generate_file_name = prefix // trim(generate_file_name) // postfix
end function

end module
