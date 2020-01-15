program inttostring
    implicit none
    character(len=1024) :: filename
    character(len=1024) :: format_string
    integer :: i
    real::r
    logical::facil=.true. 

    if (facil) then
     i=1 
     write(filename,*) i
     filename=adjustl(filename)
     write(*,*) trim(filename)
     r=2.2
     write(filename,*) r
     filename=adjustl(filename)
     write(*,*) trim(filename)
    else

    do i=1, 10
        if (i < 10) then
            format_string = "(A5,I1,F4.2,A4)"
        else
            format_string = "(A5,I2,F4.2,A4)"
        endif
        r=2.234856
        write (filename,format_string) "hello", i,r,".dat"
        write(*,*) trim(filename)
    enddo
    endif

end program inttostring
