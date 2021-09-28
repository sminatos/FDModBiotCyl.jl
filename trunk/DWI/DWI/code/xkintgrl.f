      subroutine getf(farray,gfcn,xkernal,nf,nz,nt,dk,znorm,
     *                nrecvrs)

CCCC---------------  D. R. Burns, January, 1990


c read in data from trace_green and compute wavenumber integral

      complex xkernal(nz),gfcn(10,nt),scale,k,lr,w,vfluid,temp,i0,i1
      dimension farray(nt/2),znorm(nrecvrs)
       real kreal

       common/xlayers/ r(10),vp(10),vs(10),rho(10),qp(10),qs(10)

      parameter (scale=(.63661977,0.),pi=3.1415927,twopi=6.2831853)
      data ndfft/1/,zero/0./,isign/-1/,iform/1/

 1000 format (i1)

      nf=1
      dz=1/(nz*dk)

c ----- read in the frequency and the number of k values
c ----- for this frequency

   40 read (3,end=90) farray(nf),nk

        if (nf.ne.1.and.farray(nf).le.farray(nf-1)) then
         write (6,4000) farray(nf-1),farray(nf)
 4000    format (' frequencies not increasing ',2(g9.4,1x)/)
         nf=nf-1
         irep=1
        else
         irep=0
        end if

c ----- input the terms for the integral for this f
c SM: reading xkernal==a1prme from 3 (written in tracegreen_table.f)
      read (3,end=990) (xkernal(i),i=1,nk)
        if(nf.eq.1) then
          write(64,*) farray(nf),nk,(xkernal(i),i=1,nk)
        endif


        if (irep.eq.1) go to 40

c ----- multiply the A1prime terms by -2i/pi so they are scaled properly
c ----- with the primary excitation (to be calculated later)
c ----- (Note that here we really only put in 2/pi because there is a 
c ----- factor of i out in front of the integral.)
cSM A1prime terms --> xkernal?
C     SM: So, where (-)sign goes then? (is this the source of flipped sign of spect in trace_time.f?)
C       : or isign=-1 in fourt in this code?        

   51 do 61 iscale=1,nk
      xkernal(iscale)=xkernal(iscale)*scale
   61 continue


c ----- do the integral with an fft


        if (nk.gt.nz/2.) then
         write (6,5000)
 5000    format ('** WARNING nk too big')
        end if

c ----- zero out second part of array
      do 70 i=nk+1,nz
      xkernal(i)=(0.,0.)
   70 continue

c ----- make the wavenumber data symmetric for the fft
c ----- only folding is necessary in the k dimension
      do 80 i=2,nk
      xkernal(nz+2-i)=xkernal(i)
   80 continue

c ----- do the Fourier transform
      call fourt(xkernal,nz,ndfft,isign,iform,zero)

        do 9090 ipts=1,nrecvrs
          iz=int(znorm(ipts)/dz/twopi)+1
          gfcn(ipts,nf)=xkernal(iz)*dk
9090    continue

c ----- prepare to read terms for next f
      nf=nf+1
        if (nf.gt.nt/2+1) then
        write (6,6000)
 6000   format ( '** ERROR - nf too big')
        stop01
        end if
      go to 40

90      nf = nf - 1
        fmax = farray(nf)
        return

c ----- incomplete data for some frequency
  990 write (6,99000) farray(nf),i,nk
99000 format ('** WARNING Incomplete data for frequency ',f8.4/,
     *        '   reading point ',i4,' out of ',i4)
      nk=i
      go to 51

      end
