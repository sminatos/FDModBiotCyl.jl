CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC  TRACE_TIME.F  CCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c   D. Burns                      12/6/89
c                          mods:   10/4/91
c
c   based on "ltseries.f" by K. Tubman, ERL, MIT @1984
c
c
c ----- do the k integrations for a given source-receiver separation using
c ----- ffts then use another fft to get the time series
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine trace_time

c ----- work variables


      complex ci,tshft,source(1024),spect(1024),gfcn(10,1024),cn,wo,
     *  a1prme(65536),scale,primex(1024)

      real klast

c ----- work arrays

      dimension wr(1024),yorig(1024),y(1024),fsave(256),
     *          gorig(1024),znorm(10)

c------ amplitude scaling array
       dimension ampscal(20)
c ----- fft array
      dimension nn(1)
c ----- parameter arrays

        common/xlayers/ r(10),vp(10),vs(10),rho(10),qp(10),qs(10)


      parameter (ci=(0.,1.),scale=(.63661977,0.),pi=3.1415927,
     *           twopi=6.2831853)

c ----- nt - number of time points - dimension of spect,gfcn,
c -----                              primex,wr,yorig,y2,y,time
c ----- nt/2 is the dimension of fsave
c ----- nz - number of z points - dimension of a1prme

        character*80 fout

        common/timing/fout,z,nsrcpts,dzsource,isrc,f0,alpha
	common /hayashi/t0,tlast

      data nt/1024/,nz/65536/

        write(6,*) 'params:',fout,z,nsrcpts,dzsource,isrc,f0,alpha

c ----- variables needed for the fft routine

      data ndfft/1/,iform/1/


c ----- initalize some variables
      td=0.
      ts=-99.
c
c****************************************
c************READ IN********************
c

30    continue

       write(6,*) 'in trace_time'


C  30 write (6,4000)
C4000 format (' input file name ',$)
C     read (5,*,err=30) flname

C      write(6,*) 'enter output filename'
C       read(5,*) fout
C       write(6,*) fout


C        open (3,file=flname,err=999,form='unformatted')
         open (3,form='unformatted')

          open(12,file=fout,form='formatted')


   35 rewind (3)

C  40 write (6,5000)
C5000 format (' enter source-receiver separation ',$)
C        read(5,*) z

C      write(6,*) 'enter number of receivers'
C       read(5,*) nsrcpts

C       write(6,*) 'enter receiver offsets in feet'
C        read(5,*) dzsource


C      write(6,*) 'enter radial distance (normalized) of' 
C      write(6,*) 'receiver position (0=bhole axis; 1=wall)' 
C         read(5,*) radout
          radout=0

C        write(6,*) 'displ (=1) or pressure(.ne.1)?'
C          read(5,*) idispl
           idispl=0

      read (3) rnorm,conv,wi,dk,klast,tstart,tend,fref
      write(6,*) rnorm,conv,wi,dk,klast,tstart,tend,fref
      nlayrs=1
   50 read (3) r(nlayrs),vp(nlayrs),vs(nlayrs),
     *                      rho(nlayrs),qp(nlayrs),qs(nlayrs)
      write(6,*) r(nlayrs),vp(nlayrs),vs(nlayrs),
     *                      rho(nlayrs),qp(nlayrs),qs(nlayrs)
        if (r(nlayrs).eq.0) go to 60
       nlayrs=nlayrs+1
       go to 50

   60 znorm1=z/rnorm
        dznorm=dzsource/rnorm


c      write(6,*) rnorm,conv,wi,dk,klast

                 do 62 ii=1, nsrcpts
62               ampscal(ii)=1.0

c*********************************************************
c**********************************************************


C        write(6,*) 'limit fmax? (1=y)'
C        read(5,*) lfmax
C          if(lfmax.eq.1) then        
C              write(6,*) 'enter fmax'
C                  read(5,*) fmax
C          endif


c************MAJOR LOOP FOR BENDER PT SOURCES*********


       do 6675  ipts=1,nsrcpts

c  set up offset value for each point
c  the first point is at an offset of "znorm"
c
c   subsequent points are at "znorm + (ipts-1)*dznorm"
c

        znorm(ipts)=znorm1+(ipts-1)*dznorm


c ----- read in the k values for each frequeny and do the integral
c ----- to get a spectral value for each frequency



6675       CONTINUE

C SM: getf is in xkintgrl.f  (get gfcn from a1prme: coefficient A1' in Tubman1984?)         
      call getf(fsave,gfcn,a1prme,nf,nz,nt,dk,znorm,
     *          nsrcpts)


C      write(6,*) nf
   70 freqno=float(nf)
      do 80 if=1,nf
      wr(if)=twopi*fsave(if)
C     write(6,*) if,fsave(if),wr(if)
   80 continue

1000       format(i1) 


       do 6700 ipts=1,nsrcpts

         if(ipts.eq.1) then
c
c******************SOURCE*******************************
c*********************************************************
c
c  get source information on first pass through loop
c   skip this on subsequent passes
c

C      write(26,*) dtz

      call adsource(tstart,tend,nf,wr,isrc,f0,alpha,
     +              nt,dt,dwr,source)

        write(44,*) source
C       write(26,*) dtz

c*********************************************************
c**********************FILTER*****************************


C       write(6,*) 'ifilt = ', ifilt

C       call filter(coeffs,npfilt,edge,fx,wtx,bands,maxnb,nf,
C    *        wr,nt,ifilt)


C         write(26,*) dtz
c

cc**********************************************************

           endif

c***************************
c***************************

c**********CALCULATE PRIMARY ENERGY**********
c
      do 210 izero=1,nt
      primex(izero)=(0.,0.)
  210 continue

c
c ********* loop over frequency and apply source, filter to spectra
c

      do 220 if=1,nf


      tshft=cexp(-ci*(wr(if)*(t0+td)))


c ----- generate frequency response of source (primary) excitation

      primex(if)=2./znorm(ipts)*cexp(ci*wr(if)*z/(vp(1)*(1.+1./pi/qp(1)*
     *            alog(wr(if)/twopi/fref)-ci/2./qp(1))))


c**********************************************
c**************APPLY SOURCE******************

c ----- multiply Green's function by the source spectrum


         spect(if)=gfcn(ipts,if)*source(if)*tshft
         primex(if)=primex(if)*source(if)*tshft
C        write(42,*) wr(if),gfcn(ipts,if),source(if),tshft,spect(if)
C        write(43,*) wr(if),tshft,primex(if)


c****************************************************
c*****************APPLY FILTER**********************

C       write(6,*) 'ifilt=',ifilt
C       if (ifilt.eq.1) then

C        write(6,*) 'in filter loop'
c ----- filter the signal

C       spect(if)=spect(if)*filt(if)
C       primex(if)=primex(if)*filt(if)


C       end if

  220 continue

  240   continue



CCC OPTIONAL --- OUTPUT GREEN FUNCTION, and SPECTRUM for
CCCC   PLOTTING
CCCCC

c ----- find amplitude spectrum  including source
c ----- the array y is set to this for plotting

       do 250 i=1,nf
       yorig(i)=cabs(spect(i)+primex(i))

        gorig(i)=cabs(gfcn(ipts,i))

  250 continue

       write(24,*) yorig
        write(34,*) gorig
        write(54,*) spect


c**********************************************
c***********FFT FOR TIME SERIES***************


c ----- zero out second part of array

      do 290 if=nf+1,nt
      spect(if)=(0.,0.)
      primex(if)=(0.,0.)
  290 continue

c ----- make the spectrum symmetric for the inverse fft

      do 300 if=2,nf
      spect(nt+2-if)=conjg(spect(if))
      primex(nt+2-if)=conjg(primex(if))
  300 continue


      nn(1)=nt
      isign=-1
      zero=0.
      call fourt (spect,nn,ndfft,isign,iform,zero)
      call fourt (primex,nn,ndfft,isign,iform,zero)


c****************************************************
c**********REMOVE ARTIFICIAL ATTENUATION*********


c ----- find maximum and multiply by exp(wi*t) to eliminate
c ----- artificial attenuation introduced in seis.f
c ----- also add in the primary excitation and find the 
c ----- time series

C      write(26,*) nt,nf

  310 ymax=0.


      do 320 it=1,nt

         y(it)=real(spect(it)*exp(wi*(it-1)*dt)+primex(it))*dwr
c        y(it)=real(spect(it)+primex(it))*dwr
         y(it)=y(it)*ampscal(ipts)
     
C hayashi         write(55,*) it,wi,dwr,ampscal,spect(it),primex(it)

         if (abs(y(it)).gt.ymax) ymax=abs(y(it))

  320 continue


c*************************************************
c**********SAVE FILES FOR PLOTTING**************
c
c		Trace header for 8 header file by K. Hayashi  Jan. 21 1999

		write (6,*) ipts,z+float(ipts-1)*dzsource,dzsource,t0,tlast,nt

		write (12,*) ipts
		write (12,*) int((z+float(ipts-1)*dzsource)*100.)
		write (12,*) 0
		write (12,*) 1
		write (12,*) 1
		write (12,*) 0.
		write (12,*) (tlast-t0)/float(nt)*1000.
		write (12,*) nt

          do 6696 iwrite=1,nt
C hayashi       write(12,*) iwrite,y(iwrite)
		write (12,*) y(iwrite) 
6696      continue


c         write(26,*) nt,iwrite


6700  continue

        stop

c ----- problems with input file

  999 write (6,99000) flname
99000 format('**  ERROR opening input file: ',a70)
      go to 30

      stop
      end
