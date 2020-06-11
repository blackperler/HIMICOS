
      MODULE INTEGRATION
      IMPLICIT NONE



      CONTAINS

      PURE FUNCTION trapz(x,y) ! trapz integration
      double precision trapz
      double precision, intent(in), dimension(:):: x, y
        integer:: n, i
        n = size(x)
        trapz = 0.
        do i = 2,n
         trapz = trapz + (y(i) + y(i-1))/2*(x(i)-x(i-1))
        end do
       END FUNCTION trapz
      
      PURE FUNCTION cumtrapz(x,y) ! trapz integration
        double precision, intent(in), dimension(:):: x, y
        integer:: n
        double precision, dimension(size(x)):: cumtrapz
        integer:: i
        n = size(x)
        cumtrapz(1) = 0.
        do i = 2,n
         cumtrapz(i) = cumtrapz(i-1) + (y(i) + y(i-1))/2*(x(i)-x(i-1))
        end do
       END FUNCTION cumtrapz 
      
      END MODULE INTEGRATION
