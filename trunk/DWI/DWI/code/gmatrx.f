
      subroutine gmatrx(bcond)

c      fills the G matrix

      complex einv(4,4),d(4,4),g(4,4),w,vp(10),vs(10)
      complex k,bcond(6)
      common/layer/vp,vs,rho(10),r(10)
      common/values/ w,k,nlayrs

      call dmatrx(d,nlayrs,r(nlayrs-1))


      if (vs(nlayrs-1).eq.0.) then
      call xeqy(d,g)

       else

         do 20 i=1,nlayrs-2
            if (vs(nlayrs-i).eq.0.) return
            call ematrx(einv,nlayrs-i)
            call cmatmult(einv,d,g)
            call xeqy(g,d)
   20    continue
       endif

       bcond(1) = g(1,1)
       bcond(2) = g(3,1)
       bcond(3) = g(4,1)
       bcond(4) = g(1,3)
       bcond(5) = g(3,3)
       bcond(6) = g(4,3)
       

      return
      end
