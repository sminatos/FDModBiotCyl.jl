C                      TRACE_MAIN.F
C
C
C                    copyright NER 1991
C
C
C      date:  November 24, 1991
C
C      written by:  D. Burns
C
C      based on work by K. Tubman, MIT (1983-1984)
C
C
C      subroutine calls:    trace_green.f     compute Green's function
C                           trace_time.f      compute time series
C
C      input file:    redirected (sample=trace.input)
C
C      output file:   named in input file (sample=trace.output) 
C
C     Comments by SM
C     : gfortran test in 2020/Jun
C     : trace_time.f was modified to force t0 to be zero, and skip using cexp()
C     : Because it produces unexpected phase rotation
C      

      
        character*80 fout
        dimension r1(10),vp1(10),vs1(10),rho1(10),qp(10),qs(10),
     +             r(10),rho(10)

        complex w,k,vp(10),vs(10)
          real i0iest,i0rest,i1iest,i1rest,k0iest,k0rest,k1iest,k1rest
          real  i0real,i0imag,i1real,i1imag,k0real,k0imag,k1real,k1imag
          complex z,i0,i1,k0,k1
          real klast

        common/tabstuff/ nxmax,nymax,stepr,stepi,xmin,xmax,ymin,ymax
C       common/tabstuff/ step,xmin,xmax,ymin,ymax
        common/tables/ i0real(450,450),i0imag(450,450),i1real(450,450),
     +        i1imag(450,450),k0real(450,450),k0imag(450,450),
     +        k1real(450,450),k1imag(450,450)
        common/timing/fout,sr1,nrec,dx,isrc,fsource,alpha
	common /hayashi/t0,tlast

c	ECR time
	common /timeplease/ telapsed
C
        data conv/0.01/, wi/0.1/, twopi/6.28318531/
C       data dk/0.01/, klast/10.0/
        data klast/10.0/

ccc   INITIALIZATION SECTION .... READ INPUT FILE PARAMETERS


           telapsed = 0.0

          read(5,*) nxmax,nymax,stepr,stepi,xmin,xmax,ymin,ymax
          write(6,*) nxmax,nymax,stepr,stepi,xmin,xmax,ymin,ymax
C         read(5,*) nx,ny,step,xmin,xmax,ymin,ymax

          write(6,*) 'green function already?'
          read(5,*) igrnflag
          write(6,*) igrnflag

          if(igrnflag.eq.0) then

C     read (5,*) rnorm,dk,klast,t0,tlast,f0,flast,fref
      read (5,*) t0,tlast,f0,flast,fref
      write (6,*) t0,tlast,f0,flast,fref
      read(5,*) dk
      write(6,*) dk

    8   continue

      rewind (3)

         nlayrs=1
   10    read (5,*,end=13) r1(nlayrs),vp1(nlayrs),vs1(nlayrs),
     *                     rho1(nlayrs),qp(nlayrs),qs(nlayrs)
C          if (vs1(nlayrs).eq.0) nflyr=nflyr+1
           if (r1(nlayrs).eq.0.) go to 15
         write (6,*) r1(nlayrs),vp1(nlayrs),vs1(nlayrs),
     *                     rho1(nlayrs),qp(nlayrs),qs(nlayrs)
         nlayrs=nlayrs+1
         go to 10

   13      if (r1(nlayrs).ne.0.) r1(nlayrs)=0.

   15    continue

c
c  if openhole then rnorm = bhole radius
c   else if casing and cement -- rnorm = original hole rad
c   
C          if(nlayrs.le.3) then
           rnorm=r1(1)
C          else
C          rnorm=r1(3)
C          endif


CC normalize radius and density

   30 do 32 j=1,nlayrs
      r(j)=r1(j)/rnorm
      rho(j)=rho1(j)/rho1(1)
   32 continue

      df=1./(tlast-t0)
c ----- add a little to flast to make sure the last point is calculated
      flast=flast+df/2.
      wref=twopi*fref
       write(6,*) 'ref', fref,twopi,wref



c  check dk value
c
c  dk=2pi/L  where L is the separation between imaginary sources
c   make sure that any arrivals from imag. source does not interfere
c   with the arrivals of interest (i.e. - is not within the time
c   window defined by t0, tlast)
c
902	    format (A)	
            read(5,902) fout
            write(6,902) fout
            read(5,*) sr1
            read(5,*) nrec
            read(5,*) dx
            read(5,*) isrc
            read(5,*) fsource
            read(5,*) alpha
            write(6,*) sr1,nrec,dx,isrc,fsource,alpha
    
             srdist=sr1+dx*(nrec-1)
 
c  start with assumed dk value of 0.1 (in data statement)
c  in this section check this value to be sure it's small enough
c  note: klast is assumed to be 10.0 this should be fine
c  klast should be large enough to handle 2pi*flast/cslow --- which is
c   about 2pi*5/4 == 8 for tube waves, and about 2pi*20/8 == 16 for 
c  higher shear modes, but there's not likely to be much energy out there
c  anyway....


            xsourceL=(3.14159265*2.0)/dk
            xLarrvl=(xsourceL-srdist)/vp1(nlayrs)

          write(6,*) '****xLarrvl, tlast****', xLarrvl, tlast

            if(xLarrvl.le.tlast) then

                 xtemp1=(tlast*1.25)*vp1(nlayrs)+srdist
                  dknew=(3.14159265*2.0)/xtemp1

                  write(6,*) 'WARNING:  dk value is too big, spatial wraparound 
     *               will occur'
                  write(6,*) 'change dk value in input file to be less than or 
     *                equal to  ',dknew
C      stop
                 write(6,*) '*****dk,dknew*****',dk,dknew

                  dk = dknew


            endif

C     rewind(5)

         write (3) rnorm,conv,wi,dk,klast,t0,tlast,fref

         do 20 i=1,nlayrs
         write (6,*) r1(i),vp1(i),vs1(i),rho1(i),
     *                          qp(i),qs(i)
   20    write (3) r1(i),vp1(i),vs1(i),rho1(i),
     *                          qp(i),qs(i)


c
c   read in bessl function tables.....
c

cDAN         open(81,file='../../testi0.table8bin',form='unformatted')
cDAN          open(82,file='../../testi1.table8bin',form='unformatted')
cDAN          open(83,file='../../testk0.table8bin',form='unformatted')
cDAN          open(84,file='../../testk1.table8bin',form='unformatted')

cDAN          open(81,file='testi0.table8bin',form='unformatted')
cDAN          open(82,file='testi1.table8bin',form='unformatted')
cDAN          open(83,file='testk0.table8bin',form='unformatted')
cDAN          open(84,file='testk1.table8bin',form='unformatted')

          open(81,file='tablei0.dan',form='unformatted')
          open(82,file='tablei1.dan',form='unformatted')
          open(83,file='tablek0.dan',form='unformatted')
          open(84,file='tablek1.dan',form='unformatted')
C            read(5,*) nx,ny,stepr,stepi,xmin,xmax,ymin,ymax

              do 121 iny = 1, nymax
             read(81) (i0real(inx,iny),i0imag(inx,iny),inx=1,nxmax)
119            continue
121           continue

              do 131 iny = 1, nymax
             read(82) (i1real(inx,iny),i1imag(inx,iny),inx=1,nxmax)
129            continue
131           continue


              do 141 iny = 1, nymax
             read(83) (k0real(inx,iny),k0imag(inx,iny),inx=1,nxmax)
139            continue
141           continue

              do 151 iny = 1, nymax
             read(84) (k1real(inx,iny),k1imag(inx,iny),inx=1,nxmax)
149            continue
151           continue



             write(6,*) 'nlayers = ', nlayrs

             call trace_green(vp1,vs1,qp,qs,df,f0,flast,wref,
     +                klast,dk,rnorm,nlayrs,rho,r)

        	write(*,*)'in main telapsed =',telapsed

222          call trace_time

          else


C         read (5,*) rnorm,dk,klast,t0,tlast,f0,flast,fref
          read (5,*) t0,tlast,f0,flast,fref

         nlayrs=1
1210      read (5,*,end=1213) r1(nlayrs),vp1(nlayrs),vs1(nlayrs),
     *                     rho1(nlayrs),qp(nlayrs),qs(nlayrs)
C          if (vs1(nlayrs).eq.0) nflyr=nflyr+1
           if (r1(nlayrs).eq.0.) go to 1215
         nlayrs=nlayrs+1
         go to 1210

1213      if (r1(nlayrs).ne.0.) r1(nlayrs)=0.

1215    continue

            read(5,*) fout
            read(5,*) sr1
            read(5,*) nrec
            read(5,*) dx
            read(5,*) isrc
            read(5,*) fsource
            read(5,*) alpha
            write(6,*) fout,sr1,nrec,dx,isrc,fsource,alpha

            call trace_time
c     SM: trace_time: constructs source spectrum (source), multiply with Green's function (gfcn)
c       : choose displacement or pressure (default), at axis(default) or wall     
c       : tracegreen_table.f to calc-write xkernals==A1', then xkintgrl.f load them and output gfcn
            
         endif


          stop
          end