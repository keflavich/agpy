c FILE: fplfit.f
      subroutine plfit(x,nosmall,ksa,lx)
c     internal section of plfit
c      requires that x be sorted!
      integer lx,nosmall
      real*8 x(lx)
      real*8 ksa(lx)
cf2py intent(in) :: x, nosmall
cf2py intent(hide) :: lx
cf2py intent(out) :: ksa
      integer j,n,k,flag
      real*8 lzs,cx,cf,ks,ksp,xmin,a
      xmin = 0.0
      flag = 0
c     write(*,*) "Starting length ",lx," run"
      do 100 i=1,lx-1
            if ( x(i) .lt. xmin ) write(*,*) "WARNING X WAS NOT SORTED!"
            if ( x(i) .eq. xmin ) then
                ksa(i)=0
                goto 100
            endif
            xmin = x(i)
            j=i
            lzs=0
            cx=0
            cf=0
            ks=0
            do 200 while (x(j) .ge. xmin .and. j.lt.lx) 
                n=j-i
                lzs = lzs + log(x(j)/xmin)
                if (lzs+1.eq.lzs .and. flag.eq.0) then 
                    flag=1
                    write(*,*) "Debug:",lzs,x(j),xmin,log(x(j)/xmin)
                endif
                j = j+1
200         continue            
            a = float(n) / (lzs)
            if (nosmall.gt.0) then
                if ((a-1.0)/sqrt(float(n)) .gt. 0.1) then
c                   write(*,*) "Exiting nosmall - n=",n
                    return
                endif
            endif
            ksp = 0
c           if (mod(i,100).eq.0) write(*,*) "i=",i," a=",a,"n=",
c    &       n,"lzs=",lzs
            do 300 k=0,n 
                cx = float(k) / float(n)
                cf = 1.0 - (xmin/x(k+i))**a
                ks = abs(cx-cf)
                if (ks.gt.ksp)  ksp = ks 
300         continue            
            ksa(i) = ksp
100   continue
      return
      END
c END FILE fplfit.f
