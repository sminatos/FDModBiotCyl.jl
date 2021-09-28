CCCCCCCCCCCCCCCCCCCCCCCC
CCC
CCC    SOURCE.F
CCC
CCC    computation of source function for synthetic
CCC   
CCC     D. Burns Nov, 1991
CCC              Feb, 1992
CCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine adsource(tstart,tend,maxspt,nf,w,isrc,f0,
     +                    alphsv,nt,dt,dwr,source)

      complex source(1024),ci,cn,wo,almio2
      dimension w(1024)

      parameter (ci=(0.,1.),pi=3.1415927,twopi=6.2831853)

 1000 format (i1)

c ----- isrc - type of source
c -----    0 - Tsang and Rader (1979,Geophysics,44,1706-1720)
c -----    1 - Ricker (1977,Transient waves in visco-elastic media,Elsevier)
        write(6,*) 'in source'
        write(6,*) tstart,tend,maxspt,nf,isrc,f0,alphsv


        deltat=tend-tstart

          t0=tstart
          t1=tend


          alpha=alphsv*f0*2.
          w0=f0*twopi
          dfreq=1./(t1-t0)
          dwr=twopi*dfreq
          dt=1./(nt*dfreq)
          wo=-ci*w0

             if (isrc.eq.0) then
              cn=8.*alpha*w0*wo/(alpha*alpha+2*wo*alpha)**2
             else
              cn=exp(-1.)
             end if



       do 2000 jj=1,nf

        if (isrc.eq.0) then
         source(jj)=8.*alpha*w0*(-ci*w(jj))/
     +              ((alpha-ci*w(jj)**2)+w0*w0)**2
        else
         e=w(jj)/w0
         e2=e*e
         rickr=4.40606148*e
         source(jj)=cmplx(e2*exp(-e2),0.)*cexp((0.,1.)*(rickr))
        endif


c ----- tune the source (taper out the low frequency portion)

C       if((itune.eq.1).and.(w(jj).le.w0))source(jj)=
C    *               source(jj)*cabs(source(jj))/cn


2000     continue

      return
      end
