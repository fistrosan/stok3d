MODULE CONS_PARAM
!
! J M Borrero
! Jan 7, 2007
! HAO-NCAR for HMI-Stanford
!
INTEGER,          PARAMETER  :: DP = KIND(1.D0)
INTEGER,          PARAMETER  :: SP = KIND(1.)
REAL(DP),   PARAMETER :: DPI = 3.141592653589793238462643
REAL(SP),   PARAMETER :: SPI = 3.141592653589793238462643
REAL(DP),   PARAMETER :: DD2R = DPI/180D0                                   ! deg -> rad (double precision)
REAL(SP),   PARAMETER :: SD2R = SPI/180.                                    ! deg -> rad (single precision)
REAL(DP),   PARAMETER :: LIGHT = 2.99792458E+10                             ! Speed of light (cm/s)
REAL(DP),   PARAMETER :: KBOL = 1.38056503E-16                              ! Boltzmann constant in (erg/K = cm^2 g/s^2 K)
REAL(DP),   PARAMETER :: MHYD = 1.67262158E-24                              ! Proton mass (g)
REAL(DP),   PARAMETER :: MELE = 9.10938215E-28                              ! Electron mass (g)
REAL(DP),   PARAMETER :: QELE = 4.80320425E-10                              ! Electron charge (erg*cm)^(1/2) = statCoulomb
REAL(DP),   PARAMETER :: EVOLT = 1.60217656E-12                             ! Electronvolt (ergs = cm^2 g/s^2)
REAL(DP),   PARAMETER :: MAMU = 1.66053892E-24                              ! Atomic mass unit (g): 1/12 weight of carbon atom
REAL(DP),   PARAMETER :: HPLA = 6.62606957E-27                              ! Planck's constant (erg s = cm^2 g/s)
REAL(DP),   PARAMETER :: CSAHA1 = HPLA**3D0/(2D0*DPI*MELE*KBOL)**(3D0/2D0)  ! Constant # 1 in Saha Equation
REAL(DP),   PARAMETER :: CSAHA2 = EVOLT/KBOL                                ! Constant # 2 in Saha Equation
REAL(DP),   PARAMETER :: RBOHR = HPLA**2D0/(4D0*DPI**2D0*MELE*QELE**2D0)    ! Bohr's radius (cm)
REAL(DP),   PARAMETER :: GRAV = 2.7414E4                                    ! Solar surface gravity (cm/s^2)
END MODULE CONS_PARAM
