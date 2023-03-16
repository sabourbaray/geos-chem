!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: stub_fullchem_csoa_AerosolFuncs.F90
!
! !DESCRIPTION: Stub routines for complex SOA aerosol rate-law functions etc. 
!\\
!\\
! !INTERFACE: 
!
MODULE fullchem_csoa_AerosolFuncs
!
! !USES:
!
  USE gckpp_Precision

  IMPLICIT NONE
!
! !DEFINED PARAMETERS
!
  INTEGER,  PARAMETER :: MHC   = 11 ! max # HCs
  INTEGER,  PARAMETER :: MSV   = 5  ! max # lumped semivols
  INTEGER,  PARAMETER :: MPROD = 4  ! max # volatility products
  INTEGER,  PARAMETER :: MNOX  = 3  ! max # NOx levels/oxidants
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS

 SUBROUTINE SOA_Kinetic_Rates( I, J, L, TK, KO3, KOH, KNO3, KRO2NO, KRO2HO2 )
   !
   USE gckpp_Global, ONLY : dp. MHC
   !
   INTEGER,  INTENT(IN)  :: I
   INTEGER,  INTENT(IN)  :: J
   INTEGER,  INTENT(IN)  :: L
   REAL(dp), INTENT(IN)  :: TK
   !
   REAL(dp), INTENT(OUT) :: KO3(MHC)
   REAL(dp), INTENT(OUT) :: KOH(MHC)
   REAL(dp), INTENT(OUT) :: KNO3(MHC)
   REAL(fp), INTENT(OUT) :: KRO2NO
   REAL(dp), INTENT(OUT) :: KRO2HO2
   !
 END SUBROUTINE SOA_Kinetic_Rates
 !
 SUBROUTINE SOA_Equilib_Coeffs( I, J, L, TK, KOM )
   !
   USE gckpp_Global, ONLY : dp, MPROD, MSV
   !
   INTEGER,  INTENT(IN)  :: I               ! Longitude index
   INTEGER,  INTENT(IN)  :: J               ! Latitude index
   INTEGER,  INTENT(IN)  :: L               ! Altitude index
   REAL(dp), INTENT(IN)  :: TK              ! Temperature [K]
   !
   REAL(dp), INTENT(OUT) :: KOM(MPROD,MSV)  ! Equilibrium gas-aerosol
   !
 END SUBROUTINE SOA_Equilib_Coeffs
 !
 FUNCTION SemiVolOxidation( kPhot, a_r, conc ) RESULT( k )
   !
   USE gckpp_Global, ONLY : dp, MPROD, MSV
   !
   REAL(dp), INTENT(IN) :: kPhot  ! Photo-oxidation rate      [cm3/molec/s]
   REAL(dp), INTENT(IN) :: a_r    ! Activation energy / R     [K          ]
   REAL(dp), INTENT(IN) :: conc   ! Oxidant concentration     [molec/cm3  ]
   REAL(dp)             :: k      ! Oxidation rate of semivol [1/s        ]
   !
 END FUNCTION SemiVolOxidation
 !
 FUNCTION SemiVol13TempAdjKOM( kom_ref, heat_vapor ) RESULT( kom )
   !
   ! Multiplies the reference KOM (equilibrium coefficient)
   ! by the heat of vaporization and temperature terms.
   ! For lumped semivolatiles 1 and 3.
   ! 
   USE gckpp_Global, ONLY : dp
   !
   REAL(dp), INTENT(IN) :: kom_ref
   REAL(dp), INTENT(IN) :: heat_vapor
   !
   REAL(dp)             :: kom
   !
 END FUNCTION SemiVol13TempAdjKOM
 !
 FUNCTION SemiVol45TempAdjKOM( kom_ref, heat_vapor ) RESULT( kom )
   !
   USE gckpp_Global, ONLY : dp
   !
   REAL(dp), INTENT(IN) :: kom_ref
   REAL(dp), INTENT(IN) :: heat_vapor
   !
   REAL(dp)             :: kom
   !
 END FUNCTION SemiVol45TempAdjKOM
!EOC
END MODULE fullchem_csoa_AerosolFuncs
