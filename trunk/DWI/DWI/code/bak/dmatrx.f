
      subroutine dmatrx(d,n,radius)

c ----- fills the D matrix, D(n,radius)

      complex i0,i1,k0,k1,l,m,lr,mr
      complex w,p,k,kk,ww,vs2
      complex d(4,4),vp(10),vs(10)
      common/layer/ vp,vs,rho(10),r(10)
      common/values/ w,k,nlayrs
c	ECR time
	common /timeplease/ telapsed
	real*4 tarray(2)
	static tinc
 
c ----- compute p and s wave numbers l and m.
 
      vs2=vs(n)*vs(n)
      kk=k*k
      ww=w*w
      p=rho(n)*(2.*kk*vs2-ww)
c ----- compute p and s wave numbers l and m.
      l=csqrt(kk-ww/vp(n)**2)
      lr=l*radius
        if (vs(n).ne.0.)  then
        m=csqrt(kk-ww/vs2)
        mr=m*radius
        end if
	
	ajunktime = dtime(tarray) 
      call bessl(lr,i0,i1,k0,k1,0)
	telapsed = telapsed + dtime(tarray)

 
      d(1,2)=l*i1
      d(3,2)=p*i0-2.*rho(n)*vs2*l/radius*i1

c ----- The following terms are not calculated at r=0.
c ----- They must be zero there to satisfy the boundary condition
c ----- at the center of the borehole.

        if (radius.eq.0.) return

      d(2,2)=k*i0
      d(4,2)=2.*rho(n)*vs2*k*l*i1
 
      d(1,1)=-l*k1
      d(2,1)=k*k0
      d(3,1)=p*k0+2.*rho(n)*vs2*l/radius*k1
      d(4,1)=-2.*rho(n)*vs2*k*l*k1
 
c ----- The following terms are not calculated in any fluid layer,
c ----- since B=B'=0 for fluids.
 
        if (vs(n).eq.0.) return
 
	ajunktime = dtime(tarray) 
      call bessl(mr,i0,i1,k0,k1,0)
	telapsed = telapsed + dtime(tarray)

	if(telapsed.gt.tinc)then
	tinc = tinc + 30
	write(*,*)'telapsed =',telapsed
	endif


 
      d(1,3)=-k*k1
      d(2,3)=m*k0
      d(3,3)=2.*rho(n)*vs2*k*(m*k0+k1/radius)
      d(4,3)=-p*k1
 
      d(1,4)=-k*i1
      d(2,4)=-m*i0
      d(3,4)=-2.*rho(n)*vs2*k*(m*i0-i1/radius)
      d(4,4)=-p*i1

      return
      end
