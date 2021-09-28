      subroutine analinv(d,a,det,r,w,rho)
      complex a(4,4),d(4,4),det
      complex w,k


C
C  DRB  ....... USE ANALYTIC INVERSE FROM SCHMITT AND BOUCHON (1985)
C               GEOPHYSICS (P1777)




        a(1,1) = d(3,2)
 	a(2,1) = -d(3,1)
	a(3,1) = d(3,4)
	a(4,1) = -d(3,3)
	a(1,2) = d(4,2)
	a(2,2) = -d(4,1)
	a(3,2) = d(4,4)
	a(4,2) = -d(4,3)
	a(1,3) = -d(1,2)
	a(2,3) = d(1,1)
	a(3,3) = -d(1,4)
	a(4,3) = d(1,3)
        a(1,4) = -d(2,2)
	a(2,4) = d(2,1)
	a(3,4) = -d(2,4)
	a(4,4) = d(2,3)

	do 10 i=1,4
	   do 8 j= 1,4

	     a(i,j) = a(i,j) * (r/(w**2*rho))

8          continue
10      continue

       det = w**4* rho**2/r**2

      return
      end
