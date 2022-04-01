program main
   use Shooting
   implicit none
   
   integer a, b, n
   real(8) h
   real(8), allocatable :: x(:), s(:,:), v(:,:), l(:,:), y_left(:), y_right(:)
   real(8), allocatable :: y(:)
   real(8) :: lambda, error, lambda_first
   
   a = -4
   b = 4
   n = 9

   allocate(s(n,n))
   allocate(v(n,n))
   lambda = 0

   allocate(y(n))

   !c = calc(3.0)
   h = Integration_h(a, b, n) !1,00
   call Integration_Grid(a, h, n, x) !result: -5,000 ... 5,0
   call Determine_L(s, v, n, h, l)
   call DetermineLambdaAndY(l, y, lambda)

   call Algorithm(y, h, y_left, y_right, v, lambda_first)
   print *, lambda_first !The right lambda value

end program