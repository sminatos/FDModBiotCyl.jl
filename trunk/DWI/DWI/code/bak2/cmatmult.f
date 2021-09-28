
      subroutine cmatmult(x,y,z)

      complex x(4,4),y(4,4),z(4,4)

      do 2 i=1,4
        do 1 j=1,4
         z(i,j)=x(i,1)*y(1,j)+x(i,2)*y(2,j)+x(i,3)*y(3,j)+
     *          x(i,4)*y(4,j)
1     continue
2     continue

      return
      end
