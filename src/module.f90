module Shooting
    use ExtraUsedSubroutines
    use Diagonalization

    implicit none
    save 
    contains

    real function Integration_h(a, b, n) !determine h
        integer, intent(in) :: a, b, n
        Integration_h = real((b - a) / real(n-1))
    end function

    subroutine Integration_Grid(a, h, n, x) !determine the x integration grid
        real(8), intent(in) :: h
        integer, intent(in) :: a, n
        real(8), allocatable, intent(out) :: x(:)
        integer i
        allocate(x(int(n)))

        do i = 1, int(n)
            x(i) = a + (i-1)*h
        enddo
    end subroutine

    subroutine Determine_L(s,v, n,h, l) !determine L(n)
        real(8) :: s(:,:), v(:,:), h
        real(8), allocatable, intent(out) :: l(:,:)
        integer :: n
        call S_value(s,n)
        call V_value(v,n)
        l = -(1.0 / (2.0*h**2))* s + v
    end subroutine

    subroutine DetermineLambdaAndY(l, y, lambda) !use diagonalization to determine lambda(n) and y(n)
        real(8), allocatable, intent(in) :: l(:,:)
        real(8), intent(out) :: lambda, y(:)
        real(8), allocatable:: lambda_vector(:), y_vector(:,:)
        allocate(lambda_vector(size(l(:,1))))
        allocate(y_vector(size(lambda_vector),size(lambda_vector)))
        call diagonalize(l, y_vector, lambda_vector)
        y = y_vector(:,1)
        lambda = lambda_vector(1)
    end subroutine


    subroutine Algorithm(y, h, y_left, y_right, v, lambda_first)
        real(8), intent(in) :: h, v
        real(8), allocatable, intent(out) :: y_left(:), y_right(:), y(:), lambda_first

        real(8):: y_left_accent, y_right_accent, lambda_after, error
        integer i

        error = 0.001 !tolerance level for lambda_after

        do i = 1, 100
            call ShootingMethod(y, h, y_left, y_right, v)
            call NewLambda(h, y, y_left, y_right, y_left_accent, y_right_accent, lambda_first, lambda_after)
            if ((lambda_after-lambda_first)<error)
                exit
            lambda_first = lambda_after
        enddo
    end subroutine

end module


module ExtraUsedSubroutines
    use integration_module
    implicit none
    contains

    subroutine S_value(s, n) !make S grid
        integer, intent(in) :: n
        real(8), intent(inout) :: s(:,:)
        integer :: i
        
        do i= 1, n
            s(i,i) = -2.0
            s(i,i+1)= -1.0
            s(i+1,i) = -1.0
        enddo
    end subroutine

    subroutine V_value(v, n) !make V grid
        integer, intent(in) :: n
        real(8), intent(inout) :: v(:,:)
        integer :: i

        do i= 1,n
            v(i,i) = 0.0
        enddo
    end subroutine

    subroutine NewLambda(h, y, y_left, y_right, y_left_accent, y_right_accent, lambda_first, lambda_after)
        real(8), intent(in) :: h, y_left(:), y_right(:), y(:), lambda_first
        real(8), intent(out) :: y_left_accent, y_right_accent, lambda_after

        real(8) :: v_int_left, v_int_right, difference
 
    
        call derivative(h, y_left, y_left_accent)
        call derivative(h, y_right, y_right_accent)

        call Newton_cotes(y_left, h, 1, size(y_left), v_int_left) !use trapez
        call Newton_cotes(y_right, h, 1, size(y_right), v_int_right)

        difference = 0.5 * ((y_right_accent/y_right(size(y_right)) - y_left_accent/y_left(size(y_left)) * 
        ((1/v_int_left)**2) + v_int_right/y_right(size(y_right)**2)))
        lambda_after = lambda_first - difference


    end subroutine

    subroutine derivative(h, y, y_accent)
        real(8), intent(in) :: h, y(:)
        real(8), intent(out) :: y_accent !y derivative
        y_accent = (y(size(y)) - y(size(y)-2)) /(2*h)
    end subroutine


    subroutine ShootingMethod(y, h,y_left,y_right, v)
        real(8), intent(in) :: y(:), v(:,:), h
        real(8), intent(out), allocatable :: y_left(:), y_right(:)

        integer :: number, i

        number = size(y)/2
        allocate(y_left(number))
        allocate(y_right(number))

        y_left(1) = y(1) !First 2 left and right values are known
        y_left(2) = y(2)
        y_right(1) = y(size(y))
        y_right(2) = y(size(y)-1)

        do i = 3, number !first 2 values are known, size of y_left is number. So number-2 iterations needed
            y_left(i) = -y_left(i-2) + 2*h**2 *(v(i,i) - (1/h**2)) * y_left(i+1)
            y_right(i) = -y_right(i-2) + 2*h**2 *(v(i,i)- (1/h**2))* y_right(i+1)
        enddo
    end subroutine

end module