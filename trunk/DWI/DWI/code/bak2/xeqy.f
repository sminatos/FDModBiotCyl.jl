      subroutine xeqy(x,y)


      complex x(4,4),y(4,4)

      do 2 i=1,4
         do 1 j=1,4
            y(i,j)=x(i,j)
1        continue
2     continue

      return
      end
