
          complex x,i0,i1,k0,k1
        real i0real(450,450),i0imag(450,450)
        real i1real(450,450),i1imag(450,450)
        real k0real(450,450),k0imag(450,450)
        real k1real(450,450),k1imag(450,450)

        character*80 fouti0,fouti1,foutk0,foutk1

   
             write(*,*) 'enter real argument start value (lo)'
             read(*,*) xrlo

             write(*,*) 'enter imag argument start value (lo)'
             read(*,*) xilo

             write(*,*) 'enter stepsize (real)'
             read(*,*) stepreal
             write(*,*) 'enter stepsize (imag)'
             read(*,*) stepimag

             write(*,*) 'enter nx,ny'
               read(*,*) nx,ny

        write(*,*) 'enter fout i0'
        read(*,902) fouti0
902	format (A)

        open(11,file=fouti0,form='unformatted')

        write(*,*) 'enter fout i1'
        read(*,902) fouti1

        open(12,file=fouti1,form='unformatted')

        write(*,902) 'enter fout k0'
        read(*,902) foutk0

        open(13,file=foutk0,form='unformatted')

        write(*,*) 'enter fout k1'
        read(*,902) foutk1

        open(14,file=foutk1,form='unformatted')


                do 200 j = 1,ny

                   do 100  i = 1, nx


                         ximag = xilo + (j-1)*stepimag
                         xreal = xrlo + (i-1)*stepreal

C                      write(*,*) i,j,xreal,ximag

C              if(xreal.le.1.e-06.and.xreal.ge.-1.0e-06.and.
C    +         ximag.le.1.e-06.and.ximag.ge.-1.0e-06) then
C                             xreal = 1.0e-4
C                             ximag = -1.0e-4
C                           endif

C                      write(*,*) i,j,xreal,ximag

                         x=cmplx(xreal,ximag)


		         call bfcns(x,i0,i1,k0,k1,0)


                         i0real(i,j) = real(i0)
                         i0imag(i,j) = aimag(i0)
                         i1real(i,j) = real(i1)
                         i1imag(i,j) = aimag(i1)
                         k0real(i,j) = real(k0)
                         k0imag(i,j) = aimag(k0)
                         k1real(i,j) = real(k1)
                         k1imag(i,j) = aimag(k1)

C                      write(*,*) xreal,ximag,i0real(i,j),i0imag(i,j)
C                      write(12,*) xreal,ximag,real(i1),aimag(i1)
C                      write(13,*) xreal,ximag,real(k0),aimag(k0)
C                      write(14,*) xreal,ximag,real(k1),aimag(k1)


100               continue
200             continue

         do 20 j=1,ny

             write(11) (i0real(i,j),i0imag(i,j),i=1,nx)
             write(12) (i1real(i,j),i1imag(i,j),i=1,nx)
             write(13) (k0real(i,j),k0imag(i,j),i=1,nx)
             write(14) (k1real(i,j),k1imag(i,j),i=1,nx)

20        continue

            stop
            end
