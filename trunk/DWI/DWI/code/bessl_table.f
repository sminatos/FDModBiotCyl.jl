
      subroutine bessl(z,i0,i1,k0,k1,ionly)


      complex z,z2o4,zt8,i0,i1,k0,k1,top,term,sum
      complex bf(100)
      parameter (gamma=.5772156649,pit3o2=4.712388981,
     *           pio2=1.5707963,twopi=6.2831853)
       common/tabstuff/ step,xmin,xmax,ymin,ymax
       common/tables/ i0real(200,200),i0imag(200,200),i1real(200,200),
     +       i1imag(200,200),k0real(200,200),k0imag(200,200),
     +       k1real(200,200),k1imag(200,200)

      data conv/1.e-12/


        if(real(z).gt.xmin.and.real(z).lt.xmax) then

        
            if(aimag(z).gt.ymin.and.aimag(z).lt.ymax) then


                    call lookup(z,i0,i1,k0,k1)

           endif

        else


c      calculate i0 and i1 using recurrence relationship eqn 9.6.26...
        if (abs(z).le.15.) then
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



       endif



      return
      end
c
c   BESSL FUNCTION LOOKUP TABLE ROUTINE
c
c   3/24/92   drb ------ read in arguments and compute then look up value, 
c           do a linear interpolation and compare


        subroutine lookup(z,i0,i1,k0,k1)

c   based on stepsize (in table), and xreal, ximag values, find the four 
c   points that surround the 'x' value, then .... linear interpolation 
c   (distance weighted average)



          common/tabstuff/ step,xmin,xmax,ymin,ymax
      common/tables/ i0real(200,200),i0imag(200,200),i1real(200,200),
     +       i1imag(200,200),k0real(200,200),k0imag(200,200),
     +       k1real(200,200),k1imag(200,200)

          real i0iest,i0rest,i1iest,i1rest,k0iest,k0rest,k1iest,k1rest
          real  i0real,i0imag,i1real,
     +          i1imag,k0real,k0imag,
     +          k1real,k1imag
          complex z,i0,i1,k0,k1


           nxmin = int(xmin/step)+1
           nxmax = int(xmax/step)
           nymax = abs(int(ymin/step))
           nymin = abs(int(ymax/step))+1


             xr = real(z)
             xi = aimag(z)

                          xrstep=xr/step
                          xistep=xi/step

                        ii1 = int(xrstep)+1
                        ii2 = ii1+1
                        jj2 = nymax + int(xistep) +1
                        jj1 = jj2-1


c        write(*,*) ii1,ii2,jj1,jj2

C       write(*,*) xval(ii1,jj1),xval(ii2,jj1),
C    +           xval(ii1,jj2),xval(ii2,jj2)




          delx = abs(xr) - abs(float(ii1-1)*step)

          dely = (xi) - (ymin+float(jj1-1)*step)

           dist1 = sqrt(delx**2+dely**2)
           dist2 = sqrt((step-delx)**2+dely**2)
           dist3 = sqrt(delx**2+(step-dely)**2)
           dist4 = sqrt((step-delx)**2+(step-dely)**2)

c          write(*,*) dist1, dist2, dist3, dist4

c
c  find i0 terms
c


           if(dist1.eq.0.0) then 
              i0rest = i0real(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 i0rest = i0real(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    i0rest = i0real(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       i0rest = i0real(ii2,jj2)
                     else
            sum1 = i0real(ii1,jj1)/dist1 + i0real(ii2,jj1)/dist2
     +              + i0real(ii1,jj2)/dist3 + i0real(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 i0rest = sum1/sum2

                  endif
                 endif
                endif
              endif

           if(dist1.eq.0.0) then 
              i0iest = i0imag(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 i0iest = i0imag(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    i0iest = i0imag(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       i0iest = i0imag(ii2,jj2)
                     else
            sum1 = i0imag(ii1,jj1)/dist1 + i0imag(ii2,jj1)/dist2
     +              + i0imag(ii1,jj2)/dist3 + i0imag(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 i0iest = sum1/sum2

                  endif
                 endif
                endif
              endif

c
c  find i1 terms
c


           if(dist1.eq.0.0) then 
              i1rest = i1real(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 i1rest = i1real(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    i1rest = i1real(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       i1rest = i1real(ii2,jj2)
                     else
            sum1 = i1real(ii1,jj1)/dist1 + i1real(ii2,jj1)/dist2
     +              + i1real(ii1,jj2)/dist3 + i1real(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 i1rest = sum1/sum2

                  endif
                 endif
                endif
              endif

           if(dist1.eq.0.0) then 
              i1iest = i1imag(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 i1iest = i1imag(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    i1iest = i1imag(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       i1iest = i1imag(ii2,jj2)
                     else
            sum1 = i1imag(ii1,jj1)/dist1 + i1imag(ii2,jj1)/dist2
     +              + i1imag(ii1,jj2)/dist3 + i1imag(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 i1iest = sum1/sum2

                  endif
                 endif
                endif
              endif



c
c  find k0 terms
c


           if(dist1.eq.0.0) then 
              k0rest = k0real(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 k0rest = k0real(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    k0rest = k0real(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       k0rest = k0real(ii2,jj2)
                     else
            sum1 = k0real(ii1,jj1)/dist1 + k0real(ii2,jj1)/dist2
     +              + k0real(ii1,jj2)/dist3 + k0real(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 k0rest = sum1/sum2

                  endif
                 endif
                endif
              endif

           if(dist1.eq.0.0) then 
              k0iest = k0imag(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 k0iest = k0imag(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    k0iest = k0imag(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       k0iest = k0imag(ii2,jj2)
                     else
            sum1 = k0imag(ii1,jj1)/dist1 + k0imag(ii2,jj1)/dist2
     +              + k0imag(ii1,jj2)/dist3 + k0imag(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 k0iest = sum1/sum2

                  endif
                 endif
                endif
              endif

c
c  find k1 terms
c


           if(dist1.eq.0.0) then 
              k1rest = k1real(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 k1rest = k1real(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    k1rest = k1real(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       k1rest = k1real(ii2,jj2)
                     else
            sum1 = k1real(ii1,jj1)/dist1 + k1real(ii2,jj1)/dist2
     +              + k1real(ii1,jj2)/dist3 + k1real(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 k1rest = sum1/sum2

                  endif
                 endif
                endif
              endif

           if(dist1.eq.0.0) then 
              k1iest = k1imag(ii1,jj1)
           else
              if(dist2.eq.0.0) then
                 k1iest = k1imag(ii2,jj1)
              else
                  if(dist3.eq.0.0) then
                    k1iest = k1imag(ii1,jj2)
                  else
                     if(dist4.eq.0.0) then
                       k1iest = k1imag(ii2,jj2)
                     else
            sum1 = k1imag(ii1,jj1)/dist1 + k1imag(ii2,jj1)/dist2
     +              + k1imag(ii1,jj2)/dist3 + k1imag(ii2,jj2)/dist4

            sum2 = 1/dist1 + 1/dist2 + 1/dist3 + 1/dist4

                 k1iest = sum1/sum2

                  endif
                 endif
                endif
              endif


c
c  construct final complex values of the bessl functions
c

        i0 = cmplx(i0rest,i0iest)
        i1 = cmplx(i1rest,i1iest)
        k0 = cmplx(k0rest,k0iest)
        k1 = cmplx(k1rest,k1iest)



C       write(*,*) xval(ii1,jj1),xval(ii2,jj1),
C    +           xval(ii1,jj2),xval(ii2,jj2)

C       write(*,*) xr,xi,ii1,ii2,jj1,jj2,delx,dely,sum1,sum2,xvalest

C       write(47,*) xval(ii1,jj1),xval(ii2,jj1),
C    +           xval(ii2,jj1),xval(ii2,jj2)

C       write(48,*) xr,xi,ii1,ii2,jj1,jj2,delx,dely,sum1,sum2,xvalest




        return
	end


