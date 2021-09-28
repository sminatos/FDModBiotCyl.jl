
      subroutine ematrx(einv,n)

c ----- fills the INVERSE of the E matrix, Einv(n,r(n),r(n-1))

      complex d(4,4),dinv(4,4),einv(4,4),w,vp(10),vs(10),det(2)
      complex k
      common/layer/vp,vs,rho(10),r(10)
      common/values/ w,k,nlayrs

      call dmatrx(d,n,r(n))
      call analinv(d,dinv,det,r(n),w,rho(n))


      call dmatrx(d,n,r(n-1))
      call cmatmult(d,dinv,einv)

      return
      end
