      FUNCTION ZBQLU01(DUMMY)
      DOUBLE PRECISION ZBQLU01,DUMMY
      call random_number(ZBQLU01)
      ZBQLU01 = ZBQLU01
      END
#############################################################

      FUNCTION ZBQLUAB(A,B)
*
*       Returns a random number uniformly distributed on (A,B)
*
      DOUBLE PRECISION A,B,ZBQLU01,ZBQLUAB

*
*       Even if A > B, this will work as B-A will then be -ve
*
      IF (A.NE.B) THEN
       ZBQLUAB = A + ( (B-A)*ZBQLU01(0.0D0) )
      ELSE
       ZBQLUAB = A
       WRITE(*,1)
      ENDIF
 1    FORMAT(/5X,'****WARNING**** (function ZBQLUAB) Upper and ',
     +'lower limits on uniform',/5X,'distribution are identical',/)
      END
####################################################################

      FUNCTION ZBQLEXP(MU)
*
*       Returns a random number exponentially distributed with
*       mean MU
*
      DOUBLE PRECISION MU,ZBQLEXP,ZBQLU01

      ZBQLEXP = 0.0D0

      IF (MU.LT.0.0D0) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      ZBQLEXP = -DLOG(ZBQLU01(0.0D0))*MU

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLEXP',/)

      END
******************************************************************

      FUNCTION ZEXPT01(lambda)

      DOUBLE PRECISION lambda, ZBQLU01, ZEXPT01
      DOUBLE PRECISION y, c1, c2, c3, tmp


      c1 = (EXP(lambda) - 1.0D0) / lambda
      c2 = LOG(2.0D0 / (1.0D0 + EXP(-lambda))) / lambda
      c3 = (1 - EXP(-lambda)) / lambda

      ZEXPT01 = c1 * ZBQLU01(0.0D0)
      if (ZEXPT01 .lt. 1.0D0) then
        return
      endif

      do
        ZEXPT01 = ZBQLU01(0.0D0)
        if (ZEXPT01 .lt. c2) then
          return
        endif
        y = 0.5D0 * ZBQLU01(0.0D0)
        if (y .gt. (1.0D0 - ZEXPT01)) then
          ZEXPT01 = 1.0D0 - ZEXPT01
          y = 1.0D0 - y
        endif
        if (ZEXPT01 .le. (c3 * (1.0D0 - y))) then
          return
        endif
        if ((y * c1) .le. (1.0D0 - ZEXPT01)) then
          return
        endif
        tmp = EXP(lambda * (1.0D0 - ZEXPT01)) - 1.0D0
        if ((y * c1 * lambda) .le. tmp) then
          return
        endif
      end do

      end