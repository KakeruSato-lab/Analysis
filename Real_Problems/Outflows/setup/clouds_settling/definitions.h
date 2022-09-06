#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              POTENTIAL
#define  FORCED_TURB             NO
#define  COOLING                 TABULATED
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 2
#define  USER_DEF_PARAMETERS     41

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          SELECTIVE

/* -- user-defined parameters (labels) -- */

#define  PAR_OPOW                0      // AGN wind/jet power   (erg / s)
#define  PAR_OSPD                1      // AGN wind speed (c) | jet bulk Lorentz factor
#define  PAR_OMDT                2      // AGN wind mass outflow rate (Msun / yr) | jet proper mass parameter
#define  PAR_OANG                3      // AGN wind/jet opening half-angle
#define  PAR_ORAD                4      // AGN wind/jet cross-sectional radius
#define  PAR_ODBH                5      // AGN wind/jet nozzle distance to black hole (tilt/precession axis)
#define  PAR_OSPH                6      // Size of spherical inner boudnary region.
#define  PAR_ODIR                7      // AGN wind/jet nozzle orientation (directino)
#define  PAR_OOMG                8
#define  PAR_OPHI                9
#define  PAR_OEFF                10
#define  PAR_ARAD                11
#define  PAR_AMBH                12
#define  PAR_AEFF                13
#define  PAR_AMLD                14
#define  PAR_ASNK                15
#define  PAR_HRHO                16
#define  PAR_HTMP                17
#define  PAR_HVX1                18
#define  PAR_HVX2                19
#define  PAR_HVX3                20
#define  PAR_HVRD                21
#define  PAR_HRAD                22
#define  PAR_WRHO                23
#define  PAR_WTRB                24
#define  PAR_WRAD                25
#define  PAR_WROT                26
#define  PAR_WX1L                27
#define  PAR_WX1H                28
#define  PAR_WX2L                29
#define  PAR_WX2H                30
#define  PAR_WX3L                31
#define  PAR_WX3H                32
#define  PAR_WVRD                33
#define  PAR_WVPL                34
#define  PAR_WVPP                35
#define  PAR_WVAN                36
#define  PAR_SGAV                37
#define  PAR_NCLD                38
#define  PAR_LOMX                39
#define  PAR_LCMX                40

/* [Beg] user-defined constants (do not change this line) */

#define  MU_NORM                 0.60364
#define  UNIT_DENSITY            ((CONST_amu) * (MU_NORM))
#define  UNIT_LENGTH             ((CONST_pc) * 1.e3)
#define  UNIT_VELOCITY           CONST_c

#define  ARTIFICIAL_VISC         NO
//#define  ASSIGN_VECTOR_POTENTIAL YES
#define  CHAR_LIMITING           YES
//#define  CT_EN_CORRECTION        YES
#define  ID_NZ_MAX               4
#define  INITIAL_SMOOTHING       NO
#define  INTERNAL_BOUNDARY       YES
#define  LIMITER                 MC_LIM
#define  RECONSTRUCT_4VEL        NO
#define  SHOCK_FLATTENING        MULTID
#define  PPM_ORDER               4
//#define  UPDATE_VECTOR_POTENTIAL YES
#define  WARNING_MESSAGES        NO

//#define  FORCED_TURB_ENERGY      2.e-3
//#define  FORCED_TURB_DECAY       0.5
//#define  FORCED_TURB_KMIN        (2 * CONST_PI)
//#define  FORCED_TURB_KMAX        (6 * CONST_PI)
//#define  FORCED_TURB_WEIGHT      0.3

#define  SHOW_TIME_STEPS         YES
#define  SHOW_TIMING             YES

/* [End] user-defined constants (do not change this line) */
