c  test bessl functions - generate look-up tables
c
c   call tubman routine - complex arguments
c

c  notes about the calcs:
c     interested in computation of i0,i1,k0,k1
c     with complex arguments which are the product of 
c     wavenumber*radius - where:
c        wavenumbers are:  l = csqrt(k2-w2/vp2)
c                          m = cswrt(k2-w2/vs2)
c        and radius values are distances from bhole center
c
c       note:  vp and vs are in kft/sec normalized by bhole radius
c              they are also complex (atten) as follows:
c
c               v = v/r*(1+(alog(w/wref)/pi*Q) - i/2*Q)
c
c     so:  for the range of Q's and freqs this could be a max of
c          about:
c		v = v/r*(1.1-0.05i)
c        and the velocities can range from:
c             vp ==>  3 - 25
c             vs ==>  1 - 15
c                 and r ===>  0.2 - 0.75 feet
c            k values vary from 0 - 10
c            w values vary from 0 - 20*2pi (about 128)
c
c   based on these values the wavenumber values have the following range
c
c        l==> 0 - 45      and m==> 0 - 128  (absolute values)
c
c       and the radius values will vary from  0  -  3 (normalized)
c
c
c    READ IN FILE WITH ARGUMENTS, THEN COMPUTE BESSL FUNCTIONS AND PLOT
c
c   3/24/92   drb ------ read in arguments and compute then look up value, 
c           do a linear interpolation and compare

          real i0real(400,400),i0imag(400,400),i1real(400,400),
     +               i1imag(400,400)
          real k0real(400,400),k0imag(400,400),k1real(400,400),
     +               k1imag(400,400)
          real i0rest, i0iest

          complex x,i0,i1,k0,k1
           character*80 fname,ftable

          common/tabstuff/ xmin,xmax,ymin,ymax,nxmax,nymax
   
         write(*,*) 'file (1) or manual(2)'
           read(*,*) itest

           if(itest.eq.1) then

             write(*,*) 'enter fname'
             read(*,*) fname
              open(3,file=fname)
               write(*,*) 'enter nlayrs'
                 read(*,*) nlayrs
               write(*,*) 'enter which layer (1-nlayr)'
                read(*,*) ilayr

            endif

             write(*,*) 'enter ftable'
             read(*,*) ftable
              open(2,file=ftable,form='unformatted')
             write(*,*) 'enter stepsize in table'
                read(*,*) stepsize
             write(*,*) 'enter minx,maxx, miny, maxy in table'
                read(*,*) xmin,xmax,ymin,ymax
             write(*,*) 'enter nx,ny'
                read(*,*) nxmax,nymax

C                 nxmax = int((xmax-xmin)/stepsize)
C                 nymax = abs(int((ymin-ymax)/stepsize))
               write(*,*) nxmax,nymax

                  do 17  iimag = 1,nymax
C                   do 16 ireal = 1,nxmax
                      read(2) (i0real(ireal,iimag),
     +                         i0imag(ireal,iimag),ireal=1,nxmax)
C                     write(*,*) dum1,dum2,i0real(ireal,iimag)
16                 continue
17                continue

             if(itest.eq.2) then
                 ilayr = 1
                 nlayrs = 1
              endif


25              continue

                 do 90 ii=1,nlayrs

                  if(itest.eq.1) then
                  read(3,*,end=100) dum,dum,idum,xreal,ximag
                  else
                   write(*,*) 'enter xreal,ximag'
                    read(*,*) xreal,ximag
                  endif

                    if(ii.eq.ilayr) then
                        x=cmplx(xreal,ximag)
		         call bfcns(x,i0,i1,k0,k1,0)

                     write(*,*) xreal,ximag

                         call lookup(xreal,ximag,stepsize,i0real,i0rest)
                         call lookup(xreal,ximag,stepsize,i0imag,i0iest)

                 
                 write(*,*) xreal,ximag,real(i0),aimag(i0),i0rest,i0iest
                 write(*,*) xreal,ximag,real(i1),aimag(i1),i0rest,i0iest
                 write(*,*) xreal,ximag,real(k0),aimag(k0),i0rest,i0iest
                 write(*,*) xreal,ximag,real(k1),aimag(k1),i0rest,i0iest


C                  write(21,*) xreal,ximag,real(i0),aimag(i0),i0rest,i0iest
C                  write(22,*) xreal,ximag,real(i1),aimag(i1)
C                  write(23,*) xreal,ximag,real(k0),aimag(k0)
C                  write(24,*) xreal,ximag,real(k1),aimag(k1)
                  endif
90               continue


                 go to 25

100               continue
200             continue

            stop
            end


        subroutine lookup(xr,xi,step,xval,xvalest)

c   based on stepsize (in table), and xreal, ximag values, find the four 
c   points that surround the 'x' value, then .... linear interpolation 
c   (distance weighted average)



         dimension xval(400,400), dist(4),ii2, jj2
          common/tabstuff/ xmin,xmax,ymin,ymax,nxmax,nymax
       


C          nxmin = int(xmin/step)+1
           nxmin = 1
C          nxmax = int((xmax-xmin)/step)
C          nymax = abs(int((ymin-ymax)/step))
C          nymin = abs(int(ymax/step))+1
           nymin = 1


                          xrstep=(xr-xmin)/step
                          xistep=(xi-ymax)/step

                        ii1 = int(xrstep)+1
        if(ii1.lt.nxmin) ii1=nxmin
        if(ii1.gt.nxmax) ii1=nxmax

                        ii2 = ii1+1
        if(ii2.gt.nxmax) ii2=nxmax
        if(ii2.lt.nxmin) ii2=nxmin

C                       jj2 = nymax + int(xistep) +1
                        jj2 = nymax + int(xistep)
        if(jj2.gt.nymax) jj2=nymax
        if(jj2.lt.nymin) jj2=nymin

                        jj1 = jj2-1
        if(jj1.lt.nymin) jj1=nymin
        if(jj1.gt.nymax) jj1=nymax



C       if(ii2.gt.nxmax) ii2=nxmax
C       if(ii1.lt.nxmin) ii1=nxmin
C       if(jj2.gt.nymax) jj2=nymax
C       if(jj1.lt.nymin) jj1=nymin
       
C       if(ii1.gt.nxmax) ii1=nxmax
C       if(ii2.lt.nxmin) ii2=nxmin
C       if(jj1.gt.nymax) jj1=nymax
C       if(jj2.lt.nymin) jj2=nymin

         write(*,*) ii1,ii2,jj1,jj2

        write(*,*) xval(ii1,jj1),xval(ii2,jj1),
     +           xval(ii1,jj2),xval(ii2,jj2)




          delx = abs(xr) - abs(float(ii1-1)*step+xmin)

          dely = (xi) - (ymin+float(jj1-1)*step)

           dist(1) = sqrt(delx**2+dely**2)
           dist(2) = sqrt((step-delx)**2+dely**2)
           dist(3) = sqrt(delx**2+(step-dely)**2)
           dist(4) = sqrt((step-delx)**2+(step-dely)**2)

           write(*,*) dist(1), dist(2), dist(3), dist(4)


          write(*,*) 'distweight (1) or gradient (2) or plane (3)?'
              read(*,*) imethod

          if(imethod.eq.1) then

           if(dist(1).eq.0.0) then 
              xvalest = xval(ii1,jj1)
           else
              if(dist(2).eq.0.0) then
                 xvalest = xval(ii2,jj1)
              else
                  if(dist(3).eq.0.0) then
                    xvalest = xval(ii1,jj2)
                  else
                     if(dist(4).eq.0.0) then
                       xvalest = xval(ii2,jj2)
                     else
            sum1 = xval(ii1,jj1)/dist(1) + xval(ii2,jj1)/dist(2)
     +              + xval(ii1,jj2)/dist(3) + xval(ii2,jj2)/dist(4)

            sum2 = 1/dist(1) + 1/dist(2) + 1/dist(3) + 1/dist(4)

          xvalest = sum1/sum2

                  endif
                 endif
                endif
              endif



cc
c  add section 3/16/92.....use plane through nearest three points to compute 
c     interpolated bessl function
c
c
c  first sort the dist values....find biggest dist and discard use other 3 
c
c


        else
          if(imethod.eq.2) then

         
        if(dist(1).gt.dist(2)) then
            imax = 1
          else
            imax = 2
        endif
           if(dist(imax).gt.dist(3)) then
               imax = imax
            else
               imax = 3
            endif
                if(dist(imax).gt.dist(4)) then
                  imax=imax
                 else 
                  imax = 4
                 endif

c  find plane through three points


            if(imax.eq.4) then
                  zz1 = xval(ii1,jj1)
                  zz2 = xval(ii2,jj1)
                  zz3 = xval(ii1,jj2)
                     ddx = delx
                     ddy = dely
             else
                if(imax.eq.3) then
                     zz1 = xval(ii1,jj1)
                     zz2 = xval(ii2,jj1)
                     zz3 = xval(ii2,jj2)
                     ddx = delx
                     ddy = dely
               else
                  if(imax.eq.2) then
                      zz1 = xval(ii1,jj2)
                      zz2 = xval(ii2,jj2)
                      zz3 = xval(ii1,jj1)
                     ddx = delx
                     ddy = step-dely
                  else
                       zz1 = xval(ii1,jj2)
                       zz2 = xval(ii2,jj2)
                       zz3 = xval(ii2,jj1)
                     ddx = delx
                     ddy = step-dely
                  endif
                endif
             endif

                xslope = (zz2-zz1)/step
                yslope = (zz3-zz1)/step

            write(*,*) xslope,yslope
  
             xvalest = ddx*((zz2-zz1)/step)+zz1 
                 write(*,*) ddx,xslope,delx*xslope,zz1,xvalest
             xvalest = xvalest + ddy*((zz3-zz1)/step)
                 write(*,*) ddy,yslope,dely*yslope,xvalest
                 





           else

c find normal to plane with three nearest points, find equation of plane,
c    then find z value (bessl function value) for given xr, xi input
c


        if(dist(1).gt.dist(2)) then
            imax = 1
          else
            imax = 2
        endif
           if(dist(imax).gt.dist(3)) then
               imax = imax
            else
               imax = 3
            endif
                if(dist(imax).gt.dist(4)) then
                  imax=imax
                 else 
                  imax = 4
                 endif

c  find plane through three points
c   note dx = dy = step for the plane calcs.


            if(imax.eq.4) then
                  zz1 = xval(ii1,jj1)
                  zz2 = xval(ii2,jj1)
                  zz3 = xval(ii1,jj2)
                    xx1 = float((ii1-1)*step+xmin)
                    yy1 = (ymin+float(jj1-1)*step)
                    xstep = step
                    ystep = step
             else
                if(imax.eq.3) then
                     zz2 = xval(ii1,jj1)
                     zz1 = xval(ii2,jj1)
                     zz3 = xval(ii2,jj2)
                    xx1 = float((ii2-1)*step+xmin)
                    yy1 = (ymin+float(jj1-1)*step)
                    xstep = -step
                    ystep = step
               else
                  if(imax.eq.2) then
                      zz1 = xval(ii1,jj2)
                      zz2 = xval(ii2,jj2)
                      zz3 = xval(ii1,jj1)
                    xx1 = float((ii1-1)*step+xmin)
                    yy1 = (ymin+float(jj2-1)*step)
                    xstep = step
                    ystep = -step
                  else
                       zz2 = xval(ii1,jj2)
                       zz1 = xval(ii2,jj2)
                       zz3 = xval(ii2,jj1)
                    xx1 = float((ii2-1)*step+xmin)
                    yy1 = (ymin+float(jj2-1)*step)
                    xstep = -step
                    ystep = -step
                  endif
                endif
             endif

                dz1 = (zz2-zz1)
                dz2 = (zz3-zz1)

c
c  normal equation ===== xpl*i - ypl*j + zpl*k = dvalue
c

                xpl = -dz1*ystep
                ypl = dz2*xstep
                zpl = xstep*ystep
                dvalue = xpl*xx1 - ypl*yy1 + zpl*zz1

                xvalest = (dvalue - xpl*xr + ypl*xi)/zpl

                 write(*,*) dz1,dz2,xx1,yy1,xpl,ypl,zpl,dvalue,xvalest



           endif
        endif
            
        write(*,*) xval(ii1,jj1),xval(ii2,jj1),
     +           xval(ii1,jj2),xval(ii2,jj2)

        write(*,*) xr,xi,ii1,ii2,jj1,jj2,delx,dely,sum1,sum2,xvalest

        write(47,*) xval(ii1,jj1),xval(ii2,jj1),
     +           xval(ii2,jj1),xval(ii2,jj2)

        write(48,*) xr,xi,ii1,ii2,jj1,jj2,delx,dely,sum1,sum2,xvalest

        return
	end
