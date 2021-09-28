
      subroutine bfcns(z,i0,i1,k0,k1,ionly,konly)


      complex z,z2o4,zt8,i0,i1,k0,k1,top,term,sum
      complex bf(100)
      parameter (gamma=.5772156649,pit3o2=4.712388981,
     *           pio2=1.5707963,twopi=6.2831853)
      data conv/1.e-12/

c      calculate i0 and i1 using recurrence relationship eqn 9.6.26...



C       if (abs(z).le.15.) then
        if (abs(z).le.15..and.abs(z).gt.0.) then
         n=int(8.5+1.5*abs(z))
         bf(n+2)=(0.,0.)
         bf(n+1)=(1.e-20,1.e-20)
         sum=0.
  
         do 10 i=n,1,-1
         bf(i)=float(i)*2./z*bf(i+1)+bf(i+2)
         sum=sum+bf(i)
   10    continue

c      normalize using eqn 9.6.37
         sum=exp(z)/(bf(1)+2.*(sum-bf(1)))
         i0=sum*bf(1)
         i1=sum*bf(2)

c      ... or aymptotic series eqn 9.7.1
        else
c
c  drb add if
         if(abs(z).gt.15.) then
         zt8=z*8.
         do 40 order=0,1
         sign=-1.
         amu=4.*order*order
         sum=(1.,0.)
         term=(1.,0.)
         y=1.

         do 20 x=1.,15.,2.
         term=term*(amu-x*x)/zt8/y
           if (abs(term).lt.conv) go to 30
         sum=sum+term*sign
         y=y+1.
         sign=sign*-1.
   20    continue
      
   30      if (order.eq.0) then
            i0=cexp(z)/csqrt(twopi*z)*sum
           else
            i1=cexp(z)/csqrt(twopi*z)*sum
           end if

   40 continue

c drb add this else/endif
         else
           i0 = cmplx(1.0,0.0)
           i1 = cmplx(1.0,0.0)
         endif

        end if



c      check to see if only I's are needed

      if (ionly.eq.1) return

c      calculate k0 using series...
c      ascending series eqn 9.6.13...
        if (abs(z).lt.5.) then
         k0=(log(z/2.)+gamma)*-i0
         z2o4=z*z/4.
         factk=1.
         top=(1.,0.)
         suminv=0.

         do 50 b=1.,25.
         factk=factk*b
         top=top*z2o4
         suminv=suminv+1./b
         term=top/factk/factk*suminv
         k0=k0+term
           if (abs(term).lt.conv) go to 80
   50    continue

c      ...or asymptotic series eqn 9.7.2
        else
         zt8=z*8.
         sum=(1.,0.)
         term=(1.,0.)
         y=1.

         do 60 x=1.,30.,2.
         term=-term*x*x/zt8/y
         sum=sum+term
            if (abs(term).lt.conv) go to 70
         y=y+1.
   60    continue
   70    k0=sqrt(pio2/z)/exp(z)*sum

        end if


c      find k1 using the Wronskian eqn 9.6.15
   80 k1=(1./z-i1*k0)/i0

      return
      end
