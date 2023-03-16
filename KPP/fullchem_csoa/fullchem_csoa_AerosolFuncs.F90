!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fullchem_csoa_AerosolFuncs
!
! !DESCRIPTION: Contains functions for compuing SOA reaction rates for KPP.
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
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: SOA_Kinetic_Rates
  PUBLIC :: SOA_Equilib_Coeffs
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_kinetic_rates
!
! !DESCRIPTION: Computes kinetic reaction rates for lumped semivolatile and
!  RO2 oxidation reactions.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_Kinetic_Rates( I, J, L, TK, KO3, KOH, KNO3, KRO2NO, KRO2HO2 )
!
! !USES:
!
   USE gckpp_Global
   USE gckpp_Parameters
   USE rateLawUtilFuncs, ONLY : GCARR_ac
!
! !INPUT PARAMETERS:
!
   INTEGER,  INTENT(IN)  :: I               ! Longitude index
   INTEGER,  INTENT(IN)  :: J               ! Latitude index
   INTEGER,  INTENT(IN)  :: L               ! Altitude index
   REAL(dp), INTENT(IN)  :: TK              ! Temperature [K]
!
! !OUTPUT PARAMETERS:
!
   REAL(dp), INTENT(OUT) :: KO3(MHC)  ! Semivol + O3  rates [cm3/molec/s]
   REAL(dp), INTENT(OUT) :: KOH(MHC)  ! Semivol + OH  rates [cm3/molec/s]
   REAL(dp), INTENT(OUT) :: KNO3(MHC) ! Semivol + NO3 rates [cm3/molec/s]
   REAL(dp), INTENT(OUT) :: KRO2NO    ! RO2     + NO  rate
   REAL(dp), INTENT(OUT) :: KRO2HO2   ! RO2     + HO2 rate
!
! !REMARKS:
!  References:
!  ============================================================================
!  PHOTO-OXIDATION RATE CONSTANTS OF ORGANICS come from:
!  (1 ) Atkinson, el al., Int. J. Chem.Kinet., 27: 941-955 (1995)
!  (2 ) Shu and Atkinson, JGR 100: 7275-7281 (1995)
!  (3 ) Atkinson, J. Phys. Chem. Ref. Data 26: 215-290 (1997)
!  (4 ) Some are reproduced in Table 1 of Griffin, et al., JGR 104: 3555-3567
!  (5 ) Chung and Seinfeld (2002)
!EOP
!------------------------------------------------------------------------------
!BOC
   !=========================================================================
   ! SOA_Kinetic_Rates begins here!
   !=========================================================================

   ! Rxn rate for oxidation of SV classes by O3 [cm3/molec/s]
   !                           kPhotolysis        Activ energy/R
   KO3(1)  = SemiVolOxidation(    63.668e-18_dp,  732.0_dp, C(ind_O3) )
   KO3(2)  = SemiVolOxidation(   200.000e-18_dp,  732.0_dp, C(ind_O3) )
   KO3(3)  = SemiVolOxidation(  1744.500e-18_dp,  732.0_dp, C(ind_O3) )
   KO3(4)  = SemiVolOxidation( 11650.000e-18_dp,  732.0_dp, C(ind_O3) )

   ! Rxn rate for oxidation of SV classes by OH [cm3/molec/s]
   KOH(1)  = SemiVolOxidation(    71.026e-12_dp, -400.0_dp, C(ind_OH) )
   KOH(2)  = SemiVolOxidation(   171.000e-12_dp, -400.0_dp, C(ind_OH) )
   KOH(3)  = SemiVolOxidation(   227.690e-12_dp, -400.0_dp, C(ind_OH) )
   KOH(4)  = SemiVolOxidation(   245.000e-12_dp, -400.0_dp, C(ind_OH) )

   ! Rxn rate for oxidation of SV classes by NO3 [cm3/molec/s]
   KNO3(1) = SemiVolOxidation(     6.021e-12_dp, -490.0_dp, C(ind_NO3) )
   KNO3(2) = SemiVolOxidation(    12.200e-12_dp, -490.0_dp, C(ind_NO3) )
   KNO3(3) = SemiVolOxidation(    33.913e-12_dp, -490.0_dp, C(ind_NO3) )
   KNO3(4) = SemiVolOxidation(    27.000e-12_dp, -490.0_dp, C(ind_NO3) )

   ! RO2 + NO (cf Henze et al 2008, ACP)
   KRO2NO  = GCARR_ac( 2.6e-12_dp, 350.0_dp )  ! RO2 + NO

   ! RO2 + H2O2 (cf Henze et al 2008, ACP)
   KRO2HO2 = GCARR_ac( 1.4e-12_dp, 700.0_dp )  ! RO2 + H2O2

 END SUBROUTINE SOA_Kinetic_Rates
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_equilib_coeffs
!
! !DESCRIPTION: Computes the equilibrium gas-particle coefficients for
!  lumped semivolatiles.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_Equilib_Coeffs( I, J, L, TK, KOM )
!
! !USES:
!
   USE gckpp_Global
   USE gckpp_Parameters
!
! !INPUT PARAMETERS:
!
   INTEGER,  INTENT(IN)  :: I               ! Longitude index
   INTEGER,  INTENT(IN)  :: J               ! Latitude index
   INTEGER,  INTENT(IN)  :: L               ! Altitude index
   REAL(dp), INTENT(IN)  :: TK              ! Temperature [K]
!
! !OUTPUT PARAMETERS:
!
   REAL(dp), INTENT(OUT) :: KOM(MPROD,MSV)  ! Equilibrium gas-aerosol
                                            !  partition coeff [m3/ug]
! !REMARK:
!  ACTIVATION ENERGIES come from:
!  (1) Atkinson, R. (1994) Gas-Phase Tropospheric Chemistry of Organic
!       Compounds.  J. Phys. Chem. Ref. Data, Monograph No.2, 1-216.
!  (2) They are also reproduced in Tables B.9 and B.10 of Seinfeld and
!       Pandis (1988).
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   ! Heat of vaporization (from CRC Handbook of Chemistry & Physics)
   REAL(dp), PARAMETER   :: HEAT_VAPOR  = 5.0e+3_dp

   ! Reference coefficients (only do computation once)
   REAL(dp), PARAMETER   :: KOM_REF_1_4 = 1.0_dp      / 1646.0_dp
   REAL(dp), PARAMETER   :: KOM_REF_2_4 = 1.0_dp      /   20.0_dp
   REAL(dp), PARAMETER   :: KOM_REF_1_5 = KOM_REF_1_4 *  100.0_dp
   REAL(dp), PARAMETER   :: KOM_REF_2_5 = KOM_REF_2_4 *  100.0_dp

   !=========================================================================
   ! SOA_Equilib_Coeffs begins here!
   !=========================================================================

   ! Save to variables in gckpp_Global, for the rate-law functions
   ! We will eventually pull all of these rates to KPP
   TEMP           = TK
   INV_TEMP       = 1.0_dp / TK
   TEMP_OVER_K298 = TK     / 298.0_dp
   TEMP_OVER_K300 = TK     / 300.0_dp

   !=========================================================================
   ! Equilibrium gas-particle partition coefficients of
   ! semi-volatile compounds [ug-1 m**3]
   !=========================================================================

   ! Initialize
   KOM      = 0.0_dp

   ! Semivolatile 1: MTPA, LIMO, MTPO, SESQ
   KOM(1,1) = SemiVol13TempAdjKOM(      1.0_dp, HEAT_VAPOR ) ! C* = 1
   KOM(2,1) = SemiVol13TempAdjKOM(      0.1_dp, HEAT_VAPOR ) ! C* = 10
   KOM(3,1) = SemiVol13TempAdjKOM(     0.01_dp, HEAT_VAPOR ) ! C* = 100
   KOM(4,1) = SemiVol13TempAdjKOM(     10.0_dp, HEAT_VAPOR ) ! C* = 0.1

   ! Semivolatile 2: ISOP, now skipped
   ! Semivolatile 3: BENZ,TOLU,XYLE,(NAP)
   KOM(1,3) = SemiVol13TempAdjKOM(      1.0_dp, HEAT_VAPOR ) ! C* = 1
   KOM(2,3) = SemiVol13TempAdjKOM(      0.1_dp, HEAT_VAPOR ) ! C* = 10
   KOM(3,3) = SemiVol13TempAdjKOM(     0.01_dp, HEAT_VAPOR ) ! C* = 100
   KOM(4,3) = SemiVol13TempAdjKOM(  1.0e+10_dp, HEAT_VAPOR ) ! Low NOx, nonvol

   ! Lumped semivolatile 4: POA (primary semivolatiles)
   KOM(1,4) = SemiVol45TempAdjKOM( KOM_REF_1_4, HEAT_VAPOR )
   KOM(2,4) = SemiVol45TempAdjKOM( KOM_REF_2_4, HEAT_VAPOR )

   ! Lumped semivolatile 5: OPOA (oxidized semivolatiles)
   KOM(1,5) = SemiVol45TempAdjKOM( KOM_REF_1_5, HEAT_VAPOR )
   KOM(2,5) = SemiVol45TempAdjKOM( KOM_REF_2_5, HEAT_VAPOR )

 END SUBROUTINE SOA_Equilib_Coeffs

 !###########################################################################
 !#####          Rate law functions are defined below here              #####
 !###########################################################################

 FUNCTION SemiVolOxidation( kPhot, a_r, conc ) RESULT( k )
   !
   ! Computes the reaction rate [1/s] for the oxidation of SVOA species
   ! (with oxidants O3, OH, NO3).
   !
   USE gckpp_Global, ONLY: INV_K298, INV_TEMP
   !
   REAL(dp), INTENT(IN) :: kPhot  ! Photo-oxidation rate      [cm3/molec/s]
   REAL(dp), INTENT(IN) :: a_r    ! Activation energy / R     [K          ]
   REAL(dp), INTENT(IN) :: conc   ! Oxidant concentration     [molec/cm3  ]
   !
   REAL(dp)             :: k      ! Oxidation rate of semivol [1/s        ]
   !
   k = kPhot * EXP( a_r * ( INV_K298 - INV_TEMP ) ) * conc
   !
 END FUNCTION SemiVolOxidation

 FUNCTION SemiVol13TempAdjKOM( kom_ref, heat_vapor ) RESULT( kom )
   !
   ! Multiplies the reference KOM (equilibrium coefficient)
   ! by the heat of vaporization and temperature terms.
   ! For lumped semivolatiles 1 and 3.
   !
   USE gckpp_Global, ONLY: INV_K298, INV_TEMP, TEMP_OVER_K298
   !
   REAL(dp), INTENT(IN) :: kom_ref
   REAL(dp), INTENT(IN) :: heat_vapor
   !
   REAL(dp)             :: kom
   !
   kom = kom_ref                                                             &
       * TEMP_OVER_K298                                                      &
       * EXP( heat_vapor * ( INV_TEMP - INV_K298 ) )
   !
 END FUNCTION SemiVol13TempAdjKOM

 FUNCTION SemiVol45TempAdjKOM( kom_ref, heat_vapor ) RESULT( kom )
   !
   ! Multiplies the reference KOM (equilibrium coefficient)
   ! by the heat of vaporization and temperature terms.
   ! For lumped semivolatiles 1 and 3.
   !
   USE gckpp_Global, ONLY: INV_K300, INV_TEMP, TEMP_OVER_K300
   !
   REAL(dp), INTENT(IN) :: kom_ref
   REAL(dp), INTENT(IN) :: heat_vapor
   !
   REAL(dp)             :: kom
   !
   kom = kom_ref                                                             &
       * TEMP_OVER_K300                                                      &
       * EXP( heat_vapor * ( INV_TEMP - INV_K300 ) )
   !
 END FUNCTION SemiVol45TempAdjKOM
!EOC
END MODULE fullchem_csoa_AerosolFuncs
