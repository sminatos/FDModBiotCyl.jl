CCCC
CCCC
CCCC      TRACEGREEN_TABLE.F
CCCC
CCCC       copyright NER 1992
CCCC
CCCC        D. R. Burns     April, 1992
CCCC
CCCC
C   Comments by SM


       subroutine trace_green(vp1,vs1,qp,qs,df,f0,flast,wref,
     +               klast,dk,rnorm,nlayers,rho2,rad2)


      complex bcond(6),d(4,4),vp(10),vs(10)
      complex k,w,f1,f2,temp,xkernal(65536),
     * ci,average,tcomp1,tcomp2
      real realk,klast
       dimension rho2(10),rad2(10)
      dimension vp1(10),vs1(10),rho1(10),r1(10),qp(10),qs(10)
          real i0iest,i0rest,i1iest,i1rest,k0iest,k0rest,k1iest,k1rest
          real  i0real,i0imag,i1real,
     +          i1imag,k0real,k0imag,
     +          k1real,k1imag
          complex z,i0,i1,k0,k1
       common/tabstuff/ nxmax,nymax,stepr,stepi,xmin,xmax,ymin,ymax
       common/tables/ i0real(450,450),i0imag(450,450),i1real(450,450),
     +       i1imag(450,450),k0real(450,450),k0imag(450,450),
     +       k1real(450,450),k1imag(450,450)

      common/layer/ vp,vs,rho(10),r(10)
      common/values/ w,k,nlayrs

      parameter (ci=(0.,1.),pi=3.141592654,twopi=6.283185308)

      data converge/.01/,wi/.1/
      data kernlmax/65536/

       nlayrs = nlayers

        do 11 i=1,10
          r(i) = rad2(i)
          rho(i) = rho2(i)
11      continue

       write(6,*) 'in trace_green'
       write(6,*) f0,flast,df,klast,dk


CC FREQ LOOP-------------------------------

   35 do 200 f=f0+df,flast,df
 
      ik=0
      average=(0.,0.)

        if (f.eq.0.) then
        ik=1
        xkernal(1)=(0.,0.)
        go to 110
        end if

      wr=twopi*f
      w=cmplx(wr,wi)


      do 40 j=1,nlayrs


        if (vs1(j).eq.0.) then
         vs(j)=cmplx(0.,0.)
        else
         temp=cmplx(0.,1./2./qs(j))
         vs(j)=vs1(j)/rnorm*(1.+1./pi/qs(j)*alog(wr/wref)-temp)
        end if

      temp=cmplx(0.,1./2./qp(j))
      vp(j)=vp1(j)/rnorm*(1.+1./pi/qp(j)*alog(wr/wref)-temp)
 
   40 continue

CC
CC  WAVENUMBER LOOP --------------------------------

      do 100 realk=0.,klast,dk

        if (realk.eq.0.) then
         k=cmplx(dk/10.,0.)
        else
         k=cmplx(realk,0.)
        end if

      ik=ik+1

C SM: gmatrix.f, and dmatrix.f      
        call gmatrx(bcond)
        call dmatrx(d,1,r(1))

      tcomp1=bcond(5)*bcond(3)-bcond(2)*bcond(6)
      tcomp2=bcond(4)*bcond(3)-bcond(1)*bcond(6)
      xkernal(ik)=(d(1,1)*tcomp1-d(3,1)*tcomp2)/
     *            (d(3,2)*tcomp2-d(1,2)*tcomp1)
c SM: xkernal==A1prime in Tubman1984m (2.15)?---> slightly different...

c ----- average is an average of the previous terms

        if (ik.gt.1) then
         average=(float(ik-1)*average+xkernal(ik))/float(ik)
        else
         average=xkernal(ik)
        end if


       if(cabs(xkernal(ik)/average).lt.converge.and.ik.ge.lastik)then
         go to 110
         end if

c ----- end wavenumber loop

  100 continue

CC WRITE OUTPUT
C SM: outputting xkernal==a1prme=A1'? in file
      
  110 continue
      write (3) f,ik
      write (3) (xkernal(j),j=1,ik)

      open (7,file='freqs_done')
        if (f.eq.0) rewind (7)
      write (7,2000) f,ik
 2000 format (f9.5,i6)

      lastik=ik

c ----- end frequency loop

  200 continue

      close (3)
      close (7)

      return
      end
