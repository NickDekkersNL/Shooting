module Shooting
    use ExtraUsedSubroutines
    use diagonalization

    implicit none
    save 
    contains

    real function Integration_h(a, b, n)
        integer, intent(in) :: a, b, n
        Integration_h = real((b - a) / real(n-1))
    end function

    subroutine Integration_Grid(a, h, n, x)
        real, intent(in) :: h
        integer, intent(in) :: a, n
        real, allocatable, intent(out) :: x(:)
        integer i
        allocate(x(int(n)))

        do i = 1, int(n)
            x(i) = a + (i-1)*h
        enddo
    end subroutine

    subroutine Determine_L(s,v, n,h, l)
        real :: s(:,:), v(:,:), h
        real, allocatable, intent(out) :: l(:,:)
        integer :: n
        call S_value(s,n)
        call V_value(v,n)
        l = -(1.0 / (2*h**2))* s + v
    end subroutine

    subroutine DetermineLambdaAndY(l, lambda, y)
        real, allocatable, intent(in) :: l(:,:)
        real(8), intent(out) :: lambda(:,:), y(:)
        call diagonalize(l, lambda, y)

        
    end subroutine
end module


module ExtraUsedSubroutines
    implicit none
    contains

    subroutine S_value(s, n)
        integer, intent(in) :: n
        real, intent(inout) :: s(:,:)
        integer :: i
        
        do i= 1, n
            s(i,i) = -2.0
            s(i,i+1)= -1.0
            s(i+1,i) = -1.0
        enddo
    end subroutine

    subroutine V_value(v, n)
        integer, intent(in) :: n
        real, intent(inout) :: v(:,:)
        integer :: i

        do i= 1,n
            v(i,i) = 0.0
        enddo
    end subroutine
end module