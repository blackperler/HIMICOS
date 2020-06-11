
      MODULE INTERPOLATION
      IMPLICIT NONE



      CONTAINS
      
      pure FUNCTION linspline(x,y,xx) ! linear interpolation
        double precision, dimension(:), intent(in):: x(:),y(:)
        double precision, intent(in):: xx 
        double precision:: linspline
        integer, dimension(1):: above, below
        below = maxloc(x - xx, mask = x - xx <= 0)
        if (below(1).eq.0) below = 1
        above = minloc(x - xx, mask = x - xx >= 0)
        if (above(1).eq.0) above = size(x)
        if (above(1)==below(1)) then
             linspline = y(above(1))
        else
             linspline = (y(above(1)) - y(below(1))) / &
              &  (x(above(1)) - x(below(1))) * &
              &  (xx - x(below(1))) + y(below(1))
        end if
       END FUNCTION linspline

      END MODULE INTERPOLATION

