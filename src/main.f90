program main
   use Shooting
   implicit none
   
   integer a, b, n
   real h
   real, allocatable :: x(:), s(:,:), v(:,:), l(:,:)
   real(8), allocatable :: lambda(:,:), y(:)
   a = -5
   b = 5
   n = 11

   allocate(s(n,n))
   allocate(v(n,n))
   allocate(lambda(n,n))
   allocate(y(n))

   !c = calc(3.0)
   h = Integration_h(a, b, n) !1,00
   call Integration_Grid(a,h,n,x) !-5,000 ... 5,0
   call Determine_L(s, v, n,h, l)
   call DetermineLambdaAndY(l, lambda, y)

end program