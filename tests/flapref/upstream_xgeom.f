C  Verbatim extract from XFOIL 6.99 (M. Drela, GPL; web.mit.edu/drela/Public/web/xfoil),
C  used only as the test reference for the C port of FLAP (xfoil-c/src/xfoil_ext.c).
C  Regenerate from the upstream tarball; do not edit.

      SUBROUTINE SSS(SS,S1,S2,DEL,XBF,YBF,X,XP,Y,YP,S,N,ISIDE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C----------------------------------------------------------------
C     Returns arc length points S1,S2 at flap surface break
C     locations.  S1 is on fixed airfoil part, S2 is on flap.
C     The points are defined according to two cases:
C
C
C     If DEL > 0:  Surface will be eliminated in S1 < s < S2
C
C     Returns the arc length values S1,S2 of the endpoints
C     of the airfoil surface segment which "disappears" as a
C     result of the flap deflection.  The line segments between
C     these enpoints and the flap hinge point (XBF,YBF) have
C     an included angle of DEL.  DEL is therefore the flap
C     deflection which will join up the points at S1,S2.
C     SS is an approximate arc length value near S1 and S2.
C     It is used as an initial guess for the Newton loop 
C     for S1 and S2.
C
C
C     If DEL = 0:  Surface will be created at s = S1 = S2
C
C     If DEL=0, then S1,S2 will cooincide, and will be located
C     on the airfoil surface where the segment joining the
C     point at S1,S2 and the hinge point is perpendicular to
C     the airfoil surface.  This will be the point where the
C     airfoil surface must be broken to permit a gap to open
C     as a result of the flap deflection.
C----------------------------------------------------------------
C
C---- convergence epsilon
      DATA EPS / 1.0E-5 /
C
      STOT = ABS( S(N) - S(1) )
C
      SIND = SIN(0.5*ABS(DEL))
C
      SSGN = 1.0
      IF(ISIDE.EQ.1) SSGN = -1.0
C
C---- initial guesses for S1, S2
      RSQ = (SEVAL(SS,X,XP,S,N)-XBF)**2 + (SEVAL(SS,Y,YP,S,N)-YBF)**2
      S1 = SS - (SIND*SQRT(RSQ) + EPS*STOT)*SSGN
      S2 = SS + (SIND*SQRT(RSQ) + EPS*STOT)*SSGN
C
C---- Newton iteration loop
      DO 10 ITER=1, 10
        X1  = SEVAL(S1,X,XP,S,N)
        X1P = DEVAL(S1,X,XP,S,N)
        Y1  = SEVAL(S1,Y,YP,S,N)
        Y1P = DEVAL(S1,Y,YP,S,N)
C
        X2  = SEVAL(S2,X,XP,S,N)
        X2P = DEVAL(S2,X,XP,S,N)
        Y2  = SEVAL(S2,Y,YP,S,N)
        Y2P = DEVAL(S2,Y,YP,S,N)
C
        R1SQ = (X1-XBF)**2 + (Y1-YBF)**2
        R2SQ = (X2-XBF)**2 + (Y2-YBF)**2
        R1 = SQRT(R1SQ)
        R2 = SQRT(R2SQ)
C
        RRSQ = (X1-X2)**2 + (Y1-Y2)**2
        RR = SQRT(RRSQ)
C
        IF(R1.LE.EPS*STOT .OR. R2.LE.EPS*STOT) THEN
         S1 = SS
         S2 = SS
         RETURN
        ENDIF
C
        R1_S1 = (X1P*(X1-XBF) + Y1P*(Y1-YBF))/R1
        R2_S2 = (X2P*(X2-XBF) + Y2P*(Y2-YBF))/R2
C
        IF(SIND.GT.0.01) THEN
C
         IF(RR.EQ.0.0) RETURN
C
         RR_S1 =  (X1P*(X1-X2) + Y1P*(Y1-Y2))/RR
         RR_S2 = -(X2P*(X1-X2) + Y2P*(Y1-Y2))/RR
C
C------- Residual 1: set included angle via dot product
         RS1 = ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))/RR - SIND*R1
         A11 = ((XBF-X1)*( -X1P) + (YBF-Y1)*( -Y1P))/RR
     &       + ((  -X1P)*(X2-X1) + (  -Y1P)*(Y2-Y1))/RR
     &       - ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))*RR_S1/RRSQ
     &       - SIND*R1_S1
         A12 = ((XBF-X1)*(X2P  ) + (YBF-Y1)*(Y2P  ))/RR
     &       - ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))*RR_S2/RRSQ
C
C------- Residual 2: set equal length segments
         RS2 = R1 - R2
         A21 = R1_S1
         A22 =    - R2_S2
C
        ELSE
C
C------- Residual 1: set included angle via small angle approximation
         RS1 = (R1+R2)*SIND + (S1 - S2)*SSGN
         A11 =  R1_S1 *SIND + SSGN
         A12 =  R2_S2 *SIND - SSGN
C
C------- Residual 2: set vector sum of line segments beteen the 
C-       endpoints and flap hinge to be perpendicular to airfoil surface.
         X1PP = D2VAL(S1,X,XP,S,N)
         Y1PP = D2VAL(S1,Y,YP,S,N)
         X2PP = D2VAL(S2,X,XP,S,N)
         Y2PP = D2VAL(S2,Y,YP,S,N)
C
         XTOT = X1+X2 - 2.0*XBF
         YTOT = Y1+Y2 - 2.0*YBF
C
         RS2 = XTOT*(X1P+X2P) + YTOT*(Y1P+Y2P)
         A21 =  X1P*(X1P+X2P) +  Y1P*(Y1P+Y2P) + XTOT*X1PP + YTOT*Y1PP
         A22 =  X2P*(X1P+X2P) +  Y2P*(Y1P+Y2P) + XTOT*X2PP + YTOT*Y2PP
C
        ENDIF
C
        DET =   A11*A22 - A12*A21
        DS1 = -(RS1*A22 - A12*RS2) / DET
        DS2 = -(A11*RS2 - RS1*A21) / DET
C
        DS1 = MIN( DS1 , 0.01*STOT )
        DS1 = MAX( DS1 , -.01*STOT )
        DS2 = MIN( DS2 , 0.01*STOT )
        DS2 = MAX( DS2 , -.01*STOT )
C
        S1 = S1 + DS1
        S2 = S2 + DS2
        IF(ABS(DS1)+ABS(DS2) .LT. EPS*STOT ) GO TO 11
   10 CONTINUE
      WRITE(*,*) 'SSS: failed to converge subtending angle points'
      S1 = SS
      S2 = SS
C
   11 CONTINUE
C
C---- make sure points are identical if included angle is zero.
      IF(DEL.EQ.0.0) THEN
       S1 = 0.5*(S1+S2)
       S2 = S1
      ENDIF
C
      RETURN
      END


      SUBROUTINE SCHECK(X,Y,N, STOL, LCHANGE)
C-------------------------------------------------------------
C     Removes points from an x,y spline contour wherever 
C     the size of a segment between nodes falls below a 
C     a specified threshold of the adjacent segments.  
C     The two node points defining the short segment are
C     replaced with a single node at their midpoint.
C     Note that the number of nodes may be altered by 
C     this routine.
C
C     Intended for eliminating odd "micro" panels 
C     that occur when blending a flap to a foil.
C     If LCHANGE is set on return the airfoil definition 
C     has been changed and resplining should be done.
C
C     The recommended value for STOL is 0.05 (meaning 
C     segments less than 5% of the length of either adjoining 
C     segment are removed).  4/24/01 HHY
C------------------------------------------------------
      REAL X(*), Y(*)
      LOGICAL LCHANGE
C
      LCHANGE = .FALSE.
C--- Check STOL for sanity
      IF(STOL.GT.0.3) THEN
       WRITE(*,*) 'SCHECK:  Bad value for small panels (STOL > 0.3)'
       RETURN
      ENDIF
C
 10   DO 20 I = 2, N-2
        IM1 = I-1
        IP1 = I+1
        IP2 = I+2
C
        DXM1 = X(I) - X(I-1)
        DYM1 = Y(I) - Y(I-1)
        DSM1 = SQRT(DXM1*DXM1 + DYM1*DYM1)
C
        DXP1 = X(I+1) - X(I)
        DYP1 = Y(I+1) - Y(I)
        DSP1 = SQRT(DXP1*DXP1 + DYP1*DYP1)
C
        DXP2 = X(I+2) - X(I+1)
        DYP2 = Y(I+2) - Y(I+1)
        DSP2 = SQRT(DXP2*DXP2 + DYP2*DYP2)
C
C------- Don't mess with doubled points (slope breaks)
        IF(DSP1.EQ.0.0) GO TO 20
C
        IF(DSP1.LT.STOL*DSM1 .OR. DSP1.LT.STOL*DSP2) THEN
C------- Replace node I with average of I and I+1
         X(I) = 0.5*(X(I)+X(I+1))
         Y(I) = 0.5*(Y(I)+Y(I+1))
C------- Remove node I+1
         DO L = I+1, N
           X(L) = X(L+1)
           Y(L) = Y(L+1)
         END DO         
         N = N - 1
         LCHANGE = .TRUE.
         WRITE(*,*) 'SCHECK segment removed at ',I
         GO TO 10
        ENDIF
C
 20   CONTINUE
C
      RETURN
      END ! SCHECK

