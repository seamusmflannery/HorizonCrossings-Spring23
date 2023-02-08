C-------------------------------------------------------------------------
C
C  File XSCTNS.FOR
C
C        This set of subroutines calculates the photoelectric absorption
C   cross sections for the elements H, He, C, N, O, Ne, Na, Mg, Al, Si,
C   S, Cl, A, Ca, Cr, Fe, and Ni.  The result is in cm**2/g, given the
C   photon energy in eV.  These functions are valid only over the energy
C   range 30 - 10,000 eV, but do NOT check that the input energy is
C   within the valid range.  These functions are called by TOTLX2 to
C   calculate the total effective cross section, given a set of relative
C   abundances.  They can also be used by themselves.
C
C   VERSION : 2.0
C
C   history :
C     26 may 93 : Version 2.0 - original HELIUM function replaced
C                 with updated version.  See comments above.
C
C     21 sep 93 : VAX Fortran 77 extensions removed.
C
C     22 sep 93 : bug in FANO routine (called by HELIUM) corrected
C                 that resulted in incorrect line shape and energy (19775::MBC)
C
C     23 sep 93:  updated comments (DMc)
C
C   References:
C
C      Monika Balucinska-Church and Dan McCammon
C      "Photoelectric Absorption Cross Sections with Variable Abundances"
C      Ap.J. 400, 699 (1992)
C
C      All data (except for helium) are from:
C      B. L. Henke, P. Lee, T. J. Tanaka, R. L. Shimabukuro and B. K.
C      Fujikawa, 1982, Atomic Data and Nuclear Data Tables, vol 27, p 1.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C  Real Function: ALUM
C  Source: Atomic & Nuclear Data Tables, January 1982
C
C  Description:
C      Calculates mass absorption coefficient (mu/rho) for aluminum.
C  History:  updated below L-edge (72.78 eV) - Aug 30, 1991 (MBC)
C
C  Usage:  FUNCTION ALUM(E)
C      E = Energy in eV.
C      ALUM returns mu/rho in cm**2/gm.
C
C  COMMON Blocks:
C      none
C  IMPLICIT
C      none
C  Subroutines called by ALUM:
C      none
C
C--------------------------------------------------------------------------
      FUNCTION ALUM(E)

      REAL E, Elog, X, Alum

        Elog = Alog(E)
        If(E.LT.72.78) Then
            X = 26.90487 + (3. - 9.135221)*Elog + 1.175546*Elog*Elog

        Else If(E.LT.1559.9) Then
                X = -38.1232 + 29.5161*Elog - 4.45416*Elog*Elog +
     +                        0.226204*Elog*Elog*Elog
            Else
                X = 14.6897 + 4.22743*Elog - 0.344185*Elog*Elog +
     +                        8.18542E-3*Elog*Elog*Elog
        Endif

        Alum = Exp(X)/(E*E*E)

      RETURN
      END
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C
C  REAL FUNCTION: ARGON
C  Source: Atomic & Nuclear Data Tables, January 1982
C  Description: ARGON calculates the mass absorption coefficient of argon.
C
C ****          works well for energies above 40 eV !!!
C
C  History: updated below L-edge (245 eV) - Aug 30, 1991 (MBC)
C
C  Usage: FUNCTION ARGON(E)
C      E = Energy in eV
C      ARGON returns the mass absorption cross section in cm**2/g
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C  SUBROUTINES CALLED BY ARGON:
C      NONE
C
C---------------------------------------------------------------------------
      FUNCTION ARGON(E)

      REAL E, Elog, X, Argon

        Elog = Alog(E)
        If(E.LT.245.0) Then
           X = -330.3509 + (267.7433 + 3.)*Elog - 78.90498*Elog*Elog
     &         + 10.35983*(Elog**3) - 0.5140201*(Elog**4)

            Else If(E.LT.3202.9) Then
               X = -5.71870 + (8.85812*Elog) + (-0.307357*Elog*Elog) +
     +             (0.00169351*(Elog**3)) + (-0.0138134*(Elog**4)) +
     +             (0.00120451*(Elog**5))
            Else
               X = 19.1905 + (2.74276*Elog) + (-0.164603*Elog*Elog) +
     +             (0.00165895*Elog*Elog*Elog)
        Endif

        ARGON = Exp(X)/(E*E*E)

      Return
      End
C--------------------------------------------------------------------------
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Real Function: CALC
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for calcium
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION CALC(E)
C              E = Energy in eV
C              CALC returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION CALC(E)
C
C
      REAL E,ELOG,X,CALC
C
C
      Elog = Alog(E)


      IF(E.LT.349.31) THEN
      X=-873.972 + (865.5231 + 3.)*Elog - 339.678*Elog*Elog +
     &   66.83369*(Elog**3) - 6.590398*(Elog**4) +
     &    0.2601044*(Elog**5)

            ELSE IF(E.LT.4038.1) THEN
            X=-3449.707 + (2433.409 + 3.)*Elog
     &        - 682.0668*Elog*Elog + 95.3563*(Elog**3)
     &        - 6.655018*(Elog**4) + 0.1854492*(Elog**5)

                ELSE
                X=18.89376 + (3. - 0.2903538)*Elog
     &            - 0.1377201*Elog*Elog

      ENDIF
      Calc = Exp(X)/(E*E*E)


      RETURN
      END
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C
C  REAL FUNCTION: CARBON
C  Source: Atomic & Nuclear Data Tables, Jan. 1982
C
C  Description: CARBON calculates the mass absorption cross-section of carbon.
C
C  USAGE: FUNCTION CARBON(E)
C      E = Energy in EV
C      CARBON returns the mass absorption cross-section in cm**2/g
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY CARBON:
C      NONE
C
C-------------------------------------------------------------------------
      FUNCTION CARBON(E)
C
C
      REAL E,ELOG,X,CARBON
C
C
        Elog = Alog(E)
        If(E.LT.284.0) Then
            X = 8.74161 + (7.13348*Elog) + (-1.14604*Elog*Elog) +
     +          (0.0677044*Elog*Elog*Elog)
        Else
            X = 3.81334 + (8.93626*Elog) + (-1.06905*Elog*Elog) +
     +          (0.0422195*Elog*Elog*Elog)
        Endif
        CARBON = EXP(X)/(E*E*E)
      Return
      End
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C     Real Function: CHLORN
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for chlorine
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION CHLORN(E)
C              E = Energy in eV
C              CHLORN returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION CHLORN(E)
C
C
      REAL E,ELOG,X,CHLORN
C
C
      Elog = Alog(E)


      IF(E.LT.202.0) THEN
      X=6253.247 +(3. - 8225.248)*Elog + 4491.675*Elog*Elog
     &  - 1302.145*(Elog**3) + 211.4881*(Elog**4)
     &  - 18.25547*(Elog**5) + 0.6545154*(Elog**6)

            ELSE IF(E.LT.2819.6) THEN
            X=-233.0502 + (143.9776 + 3.)*Elog
     &        - 31.12463*Elog*Elog +
     &        2.938618*(Elog**3) - 0.104096*(Elog**4)

                ELSE
                X=-23.74675 + (14.50997 + 3.)*Elog
     &            - 1.857953*Elog*Elog + 6.6208832E-2*(Elog**3)

      ENDIF
      Chlorn = Exp(X)/(E*E*E)


      RETURN
      END
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C     Real Function: CHROM
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for chromium
C
C     History: original Aug 5, 1991 (MBC)
C
C     Usage:   FUNCTION CHROM(E)
C              E = Energy in eV
C              CHROM returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION CHROM(E)
C
C
      REAL E,ELOG,X,CHROM
C
C
      Elog = Alog(E)


      IF(E.LT.598.0) THEN
      X=-0.4919405 + (12.66939 + 3.)*Elog - 5.199775*Elog*Elog +
     &  1.086566*(Elog**3) - 0.1196001*(Elog**4) +
     &  5.2152011E-3*(Elog**5)

            ELSE IF(E.LT.691.0) THEN
            X=27.29282 +(3. - 2.703336)*Elog

                ELSE IF(E.LT.5988.8) THEN
                X=-15.2525 + (13.23729 + 3.)*Elog
     &            - 1.966778*Elog*Elog + 8.062207E-2*(Elog**3)

                     ELSE
                     X=8.307041 + (2.008987 + 3.)*Elog
     &                - 0.2580816*Elog*Elog

      ENDIF
      Chrom = Exp(X)/(E*E*E)


      RETURN
      END
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C
C     Real Funcion : HELIUM
C     Source : Marr, G. V., and West, J. B., Atomic and Nuclear Data Tables,
C                (1976) 18, 497.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C
C     Description :
C     calculates mass absorption coefficient (mu/rho) in cm2/g for neutral
C     helium for the given energy in eV.
C     Cross sections come from experimental data compiled by Marr and
C     West (Atomic Data and Nuclear Data Tables (1976) 18, 497).
C     The four strongest autoionization resonances are taken into account;
C     numbers come from Oza (Phys Rev A (1986), 33, 824), and Fernley et al.
C     (J. Phys B (1987) 20, 6457).
C
C     Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV

C     Bugs :
C     if any are found please report to the authors
C
C     History :
C     this subroutine replaces the previous version of HELIUM which
C     calculated mass absoprtion coefficients based on Henke's data
C     (Henke, B. L., et al., (1982), Atomic and Nuclear Data Tables, 27, 1).
C     This version of HELIUM returns mass  absorption coefficients which
C     are in better agreement with the best experiments as well as
C     theoretical models (see Chen, W. F., Cooper, G., and Brion, C. E.,
C     (1991), Phys. Rev. A, 44, 186).  This fortran-77 version of the
C     subroutine is based on Pat Jelinsky's program written in C
C     (obtained from EUVE Archive)
C
C     History :
C     04 jan 93 : original (19775::MBC)
C
C     23 feb 93 : comments added and modified to remove VAX
C                    fortran 77 extensions
C
C     21 sep 93 : further remnants of VAX fortran 77 extensions
C                    have been removed (19775::MBC)
C
C     23 sep 93 : bug in the FANO routine has been removed (19775::MBC)
C
C     Usage : FUNCTION HELIUM(E)
C            E = Energy in eV
C
C     Common Blocks :
C           none
C
C     Implicit :
C           none
C
C     Functions called by HELIUM
C           FANO
C
C------------------------------------------------------------------------------
      FUNCTION HELIUM(E)

C    Type definitions :
C     IMPLICIT NONE
C    Global variables :
C    Structure definitions :
C    Function declarations :
      REAL FANO
C    Local constants :
      INTEGER IP
C         ( index through loop )
      PARAMETER (IP=8)
      INTEGER IF
C         ( index through loop )
      PARAMETER (IF=4)
      REAL AV
C         ( Avogadro's number )
      PARAMETER (AV=6.022045E23)
      REAL AW
C         ( atomic weight of hydrogen )
      PARAMETER (AW=4.0026E0)
C    Local variables :
       REAL LAMBDA
C          ( wavelength in Angstroms)
       REAL X
       REAL Y
       REAL SIGMA
C          ( cross section in cm2/atom)
       INTEGER I
C          ( index trough loop)
C     Import :
       REAL E
C          ( energy in eV)
C     Export :
       REAL HELIUM
C          ( cross section in cm**2/g)
C    Local data :
       REAL C1(IP)
       REAL C2(IP)
       REAL Q(IF)
       REAL NU(IF)
       REAL GAMMA(IF)

C          ( polynomial coefficients for Marr and West data)
       DATA C1 /-2.953607E1, 7.083061E0, 8.678646E-1, -1.221932E0,
     +         4.052997E-2, 1.317109E-1, -3.265795E-2, 2.500933E-3/

C          ( polynomial coefficients for Marr and West data )
       DATA C2 /-2.465188E1, 4.354679E0, -3.553024E0, 5.573040E0,
     +          -5.872938E0, 3.720797E0, -1.226919E0, 1.576657E-1/

C          ( parameters Q for resonances (Fernley et al. 1987) )
       DATA Q /2.81E0, 2.51E0, 2.45E0, 2.44E0/

C          ( parameters NU for resonances (Oza 1986) )
       DATA NU /1.610E0, 2.795E0, 3.817E0, 4.824E0/

C          ( parameters GAMMA for resonances (Oza 1986) )
       DATA GAMMA /2.64061E-3, 6.20116E-4, 2.56061E-4,
     +             1.320159E-4/


C     Start :

       LAMBDA=12398.54E0/E
       X=ALOG10(LAMBDA)

       IF(LAMBDA.GT.503.97E0) THEN
         HELIUM=0.E0

       ELSEIF(LAMBDA.LT.46.E0) THEN
         Y=0.E0
         DO 1 I=1,IP
 1         Y=Y + C2(I)*(X**(I-1))

       ELSE
         Y=0.E0
         DO 2 I=1,IP
 2         Y=Y + C1(I)*(X**(I-1))

         DO 3 I=1,IF
 3          Y=Y + ALOG10(FANO(Q(I),NU(I),GAMMA(I),LAMBDA))

       ENDIF

       SIGMA=10.E0**Y
       HELIUM=SIGMA*AV/AW

       END
C
C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C    Real Function FANO
C
C    Source :
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C
C    Description :
C     returns a Fano line profile for use by HELIUM. The form of
C     this line shape is taken from Fernley, Taylor and Seaton (above).
C
C    Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV
C
C    Bugs :
C    if any are found please report to the authors
C
C
C    History :
C     04 jan 93 : original, based on Jelinsky's program written in
C                 C (EUVE Archive) - (19775::MBC)
C
C     23 feb 93 : comments added and modified to remove VAX
C                fortran 77 extensions (19775::MC)
C
C     21 sep 93 : further remnants of VAX fortran 77 extensions
C                  removed (19775::MBC)
C
C     22 sep 93 : bug fixed in calculations of EPSI - devided by 1.807317
C                 not added (19775::MBC)
C
C
C    Common Blocks :
C           none
C    Implicit :
C           none
C
C    Functions called by FANO
C           none
C
C----------------------------------------------------------------------------
C
       FUNCTION FANO(A,B,C,LAMBDA)

C    Type definitions :
C     IMPLICIT NONE
C    Global variables :
C    Structure definitions :
C    Function declarations :
C    Local constants :
C    Local variables :
       REAL EPS
C          ( energy in Rydbergs )
       REAL EPSI
       REAL X
C          ( log_10 of wavelength in Angstroms )
C     Import :
       REAL A
C          ( Q coefficient (Fernley et al. 1987) )
       REAL B
C          ( NU coefficient (Oza 1986) )
       REAL C
C          ( GAMMA coefficient (Oza 1986) )
       REAL LAMBDA
C          ( wavelength in Angstroms )
C     Export :
       REAL FANO
C    Start :

       EPS=911.2671E0/LAMBDA
       EPSI=3.0E0 - 1.E0/(B*B) + 1.807317
       X=2.0*(EPS - EPSI)/C
       FANO=(X-A)*(X-A)/(1.0E0 + X*X)

       END
C-----------------------------------------------------------------------
C       REAL FUNCTION HYDRO(E)
C
C       DATE:3/6/84
C       AUTHOR: ANGIE BETKER
C       SOURCE: ATOMIC AND NUCLEAR DATA TABLES, JANUARY 1982
C
C        History: modified: 6/5/87 - J. Bloch - Created F77 Vax/VMS version.
C                 updated Aug 30, 1991 (MBC)
C
C
C       DESCRIPTION: HYDRO CALCULATES MU/RHO FOR HYDROGEN IN CM**2/GM
C
C       COMMON BLOCKS:
C               NONE
C
C       SUBROUTINES CALLED:
C               NONE
C
C       IMPLICIT
C               NONE
C
C---------------------------------------------------------------------------
C
        FUNCTION HYDRO(E)
C
C
        REAL E, Elog, X, HYDRO
C
C
C        IF (E.LT.180) HYDRO = (10 ** 10.11) * (E ** -3.03)
C        IF ((E.GE.180).AND.(E.LT.750)) HYDRO=(10**10.54)*(E**-3.23)
C        IF ((E.GE.750).AND.(E.LT.6000)) HYDRO=(10**10.94)*(E**-3.37)
C        IF (E.GE.6000) HYDRO = (10 ** 10.42) * (E ** -3.24)

        Elog=Alog(E)

        X=21.46941 + (3. - 2.060152)*Elog - 0.1492932*Elog*Elog
     &    + 5.4634294e-3*(Elog**3)

        HYDRO=Exp(X)/(E*E*E)

        RETURN
        END
C
C----------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Real Function: IRON
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for iron
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION IRON(E)
C              E = Energy in eV
C              IRON returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT:
C            none
C
C-------------------------------------------------------------------------
C
      FUNCTION IRON(E)
C
C
      REAL E,ELOG,X,IRON
C
C
      Elog = Alog(E)

      IF(E.LT.707.4) THEN
      X=-15.07332 + (18.94335 + 3.)*Elog - 4.862457*Elog*Elog +
     &    0.5573765*(Elog**3) - 3.0065542E-2*(Elog**4) +
     &    4.9834867E-4*(Elog**5)

          ELSE IF(E.LT.7111.2) THEN
          X=-253.0979 + (135.4238 + 3.)*Elog - 25.47119*Elog*Elog +
     &      2.08867*(Elog**3) - 6.4264648E-2*(Elog**4)

             ELSE
             X=-1.037655 + (4.022304 + 3.)*Elog
     &       - 0.3638919*Elog*Elog

      ENDIF
      Iron = Exp(X)/(E*E*E)


      RETURN
      END
C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
C     Real Function: MAGNES
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for magnesium
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION MAGNES(E)
C              E = Energy in eV
C              MAGNES returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT:
C            none
C-------------------------------------------------------------------------
C
      FUNCTION MAGNES(E)
C
C
      REAL E,ELOG,X,MAGNES
C
C
      Elog = Alog(E)


      IF(E.LT.49.45) THEN
      X=7.107172 + (0.7359418 + 3.)*Elog

            ELSE IF(E.LT.1303.4) THEN
            X=-81.32915 + (62.2775 + 3.)*Elog
     &         - 15.00826*Elog*Elog +
     &         1.558686*(Elog**3) - 6.1339621E-2*(Elog**4)

                ELSE
                X=-9.161526 + (10.07448 + 3.)*Elog
     &            - 1.435878*Elog*Elog + 5.2728362E-2*(Elog**3)

      ENDIF
      Magnes = Exp(X)/(E*E*E)


      RETURN
      END
C
C------------------------------------------------------------------------
C----------------------------------------------------------------------
C  REAL FUNCTION: NEON
C  Source: Atomic and Nuclear Data Tables, Jan. 1982
C
C  Description:  NEON calculates the mass absorption coefficient
C      for neon gas
C
C  Usage:  REAL FUNCTION NEON(E)
C      E = Energy in eV
C      NEON returns the mass absorption cross section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY NEON:
C      NONE
C
C------------------------------------------------------------------------------
      FUNCTION NEON(E)
C
C
      REAL E,ELOG,X,NEON
C
C
        Elog = Alog(E)
        If(E.LT.867.) Then
            X = -3.04041 + (13.0071*Elog) + (-1.93205*Elog*Elog) +
     +          (0.0977639*Elog*Elog*Elog)
        Else
            X = 17.6007 + (3.29278*Elog) + (-0.263065*Elog*Elog) +
     +          (5.68290E-3*Elog*Elog*Elog)
        Endif
        NEON = Exp(X)/(E*E*E)
       RETURN
      END
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C
C  REAL FUNCTION: NICKEL
C  Source: Atomic & Nuclear Data Tables, January 1982
C
C  Description: NICKEL calculates the mass absorption coefficient for nickel.
C  History: updated below L-edge (853.6 eV) - Aug 30, 1991 (MBC)
C
C  Usage: REAL FUNCTION NICKEL(E)
C      E = Energy in eV
C      NICKEL returns the mass absorption cross section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY NICKEL:
C      NONE
C
C
C---------------------------------------------------------------------------
      FUNCTION NICKEL(E)

      REAL E, ELOG, X, NICKEL

        Elog = Alog(E)

        If(E.LT.853.6) Then
        X = -7.919931 + (11.06475 + 3.)*Elog - 1.935318*Elog*Elog
     &      + 9.3929626e-2*(Elog**3)

            Else if(E.LT.8331.6) Then
            X = 3.71129 + (8.45098*Elog) + (-0.896656*Elog*Elog) +
     +      (0.0324889*Elog*Elog*Elog)
                Else
                X = 28.4989 + (0.485797*Elog)
        Endif

        NICKEL = Exp(X)/(E*E*E)

      Return
      End
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C
C  REAL FUNCTION: NITRO
C  Source: Atomic & Nuclear Data Tables, January 1982
C
C  Description:  NITRO calculates the mass absorption cross-section of
C        nitrogen as a function of energy.
C
C  Usage: REAL FUNCTION NITRO(E)
C      E = Energy in eV
C      NITRO returns the mass absorption cross-section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY NITRO:
C      NONE
C
C--------------------------------------------------------------------------
      FUNCTION NITRO(E)
C
C
      REAL E,ELOG,X,NITRO
C
C
        Elog = Alog(E)
        If(E.LT.401.) Then
            X = 9.24058 + (7.02985*Elog) + (-1.08849*Elog*Elog) +
     +      (0.0611007*Elog*Elog*Elog)
        Else
            X = -13.0353 + (15.4851*Elog) + (-1.89502*Elog*Elog) +
     +          (0.0769412*Elog*Elog*Elog)
        Endif
        NITRO = Exp(X)/(E*E*E)
      RETURN
      END
C
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C  REAL FUNCTION: OXYGEN
C
C  Source:        Atomic & Nuclear Data Tables, January 1982
C
C  Description:  OXY Calculates the mass absorption cross-section of oxygen
C      as a function of energy.
C  History: updated above K-edge (531.7 eV) - Aug 30, 1991 (MBC)
C
C  Usage: REAL FUNCTION OXYGEN(E)
C      E = Energy in eV
C      OXYGEN returns the mass absorption cross-section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C  SUBROUTINES CALLED BY OXYGEN:
C      NONE
C
C------------------------------------------------------------------------
      FUNCTION OXYGEN(E)

      REAL E, Elog, X, Oxygen

        Elog = Alog(E)
        If(E.LT.531.7) Then
            X = 2.57264 + (10.9321*Elog) + (-1.79383*Elog*Elog) +
     +          (0.102619*Elog*Elog*Elog)
        Else
            X = 16.53869 + (0.6428144 + 3.)*Elog - 0.3177744*Elog*Elog
     &          + 7.9471897e-3*(Elog**3)

        Endif

        OXYGEN = Exp(X)/(E*E*E)

      RETURN
      END
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C       FUNCTION SILICN
C
C       Source: Atomic and Nuclear Data Tables, January 1982
C
C       Description: SILICN calculates the mass absorption cross section
C                for silicon in cm**2/g.
C       History: updated Aug 30, 1991 (MBC)
C                updated March 4, 1992 (MBC)
C
C       COMMON BLOCKS:
C               NONE
C
C       IMPLICIT
C               none
C
C       SUBROUTINES CALLED:
C               NONE
C
C---------------------------------------------------------------------------
        FUNCTION SILICN(E)

          REAL E,Elog,X,Silicn
C
        Elog = Alog(E)

        If(E.LT.100.6) Then
            X = -3.066295 + (7.006248 + 3.)*Elog - 0.9627411*Elog*Elog
        Else if(E.LT.1840.0) Then
                X = -182.7217 + (125.061 + 3.)*Elog - 29.47269*Elog*Elog
     &              + 3.03284*(Elog**3) - 0.1173096*(Elog**4)
            Else
                X = -33.39074 + (18.42992 + 3.)*Elog
     &              - 2.385117*Elog*Elog + 8.887583e-2*(Elog**3)

        Endif

        SILICN = Exp(X)/(E*E*E)
        Return
        End
C-------------------------------------------------------------------------
C------------------------------------------------------------------------
C     Real Function: SODIUM
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for sodium
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION SODIUM(E)
C              E = Energy in eV
C              SODIUM returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT:
C            none
C-------------------------------------------------------------------------
C
      FUNCTION SODIUM(E)
C
C     IMPLICIT NONE
C
      REAL E,ELOG,X,SODIUM
C
C
      Elog = Alog(E)

      IF(E.LT.1071.7) THEN
      X= -2737.598 + (2798.704 + 3.)*Elog - 1009.892*Elog*Elog +
     &  87.16455*(Elog**3) + 43.20644*(Elog**4) - 15.27259*(Elog**5)
     & + 2.180531*(Elog**6) - 0.1526546*(Elog**7) +
     &   4.3137977E-3*(Elog**8)

            ELSE
            X=1.534019 + (6.261744 + 3.)*Elog
     &        - 0.9914126*Elog*Elog +
     &        3.5278253E-2*(Elog**3)


      ENDIF
      Sodium = Exp(X)/(E*E*E)


      RETURN
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     Real Function: SULFUR
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for sulfur
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION SULFUR(E)
C              E = Energy in eV
C              SULFUR returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION SULFUR(E)
C
C
      REAL E,ELOG,X,SULFUR
C
C
      Elog = Alog(E)


      IF(E.LT.165.0) THEN
      X=598.2911 + (3. - 678.2265)*Elog + 308.1133*Elog*Elog
     &  - 68.99324*(Elog**3) + 7.62458*(Elog**4)
     &  - 0.3335031*(Elog**5)

            ELSE IF(E.LT.2470.5) THEN
            X=3994.831 + (3. - 3693.886)*Elog +
     &        1417.287*Elog*Elog - 287.9909*(Elog**3) +
     &        32.70061*(Elog**4) - 1.968987*(Elog**5) +
     &        4.9149349E-2*(Elog**6)

                ELSE
                X=-22.49628 + (14.24599+ 3.)*Elog -
     &             1.848444*Elog*Elog + 6.6506132E-2*(Elog**3)

      ENDIF
      Sulfur = Exp(X)/(E*E*E)


      RETURN
      END
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------