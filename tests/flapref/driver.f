C  Test driver for the upstream FLAP (xfoil-c tests only): sets the buffer airfoil in the
C  XFOIL.INC commons, splines it, calls FLAP with RINPUT = (XH, YH, DDEF) and returns the
C  deflected buffer and the GEOPAR results.
      SUBROUTINE REFFLAP(NIN, XIN, YIN, XH, YH, DDEF, LGSYMIN,
     &                   NOUT, XOUT, YOUT, SOUT, XPOUT, YPOUT, PARS)
      INCLUDE 'XFOIL.INC'
      DIMENSION XIN(NIN), YIN(NIN), XOUT(*), YOUT(*), SOUT(*)
      DIMENSION XPOUT(*), YPOUT(*), PARS(*), RINPUT(3)
      LOGICAL LGSYMIN
C
      NB = NIN
      DO I = 1, NB
        XB(I) = XIN(I)
        YB(I) = YIN(I)
      ENDDO
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      PI = 4.0*ATAN(1.0)
      LGSYM = LGSYMIN
      LGSAME = .TRUE.
      LBFLAP = .FALSE.
      LPLOT = .TRUE.
      CH = 0.01
      XSF = 1.0
      YSF = 1.0
      XOFF = 0.0
      YOFF = 0.0
C
      RINPUT(1) = XH
      RINPUT(2) = YH
      RINPUT(3) = DDEF
      CALL FLAP(RINPUT,3)
C
      NOUT = NB
      DO I = 1, NB
        XOUT(I) = XB(I)
        YOUT(I) = YB(I)
        SOUT(I) = SB(I)
        XPOUT(I) = XBP(I)
        YPOUT(I) = YBP(I)
      ENDDO
      PARS(1) = SBLE
      PARS(2) = CHORDB
      PARS(3) = AREAB
      PARS(4) = RADBLE
      PARS(5) = ANGBTE
      PARS(6) = THICKB
      PARS(7) = CAMBRB
      PARS(8) = XBF
      PARS(9) = YBF
      PARS(10) = 0.0
      IF(LBFLAP) PARS(10) = 1.0
      PARS(11) = 0.0
      IF(LGSAME) PARS(11) = 1.0
      PARS(12) = 0.0
      IF(LGSYM) PARS(12) = 1.0
      RETURN
      END
