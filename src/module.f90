module whatever
    implicit none
    contains
    real function calc(x)
    real, intent(in) :: x
    calc = x +1

    end function



end module

