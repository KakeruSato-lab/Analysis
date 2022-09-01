/* -- Own definitions choices -- */
//TODO: Change some of the below to ini files, use https://github.com/benhoyt/inih as ini parser.

#define NOZZLE                            NOZZLE_UFO
#define NOZZLE_FILL                       NF_PRIMITIVE
#define NOZZLE_CAP                        YES

#define ACCRETION                         NO
#define SINK_METHOD                       SINK_FREEFLOW

#define FEEDBACK_CYCLE                    NO
#define FBC_DEBOOST                       NO
#define FBC_DEBOOST_MODE                  FBC_DEBOOST_MODE_3
#define FBC_NOZZLE_DT                     1.e-5

#define GRAV_POTENTIAL                    GRAV_DOUBLE_ISOTHERMAL

#define CLOUDS                            YES
#define CLOUD_REPEAT	                  NO
#define CLOUDS_MULTI   	                  NO

#define CLOUD_DENSITY                     CD_KEPLERIAN
#define CLOUD_VELOCITY                    CV_KEPLERIAN
#define CLOUD_SCALE                       CS_VELOCITY_DISPERSION
#define CLOUD_EXTRACT_ELLIPSOID           YES
#define CLOUD_EXTRACT_CENTRAL_BUFFER      YES

#define ACCRETION_OUTPUT                  NO
#define ACCRETION_OUTPUT_RATE             1.0
#define BONDI_ACCRETION_OUTPUT            NO
#define TURBULENT_BONDI_ACCRETION         NO
#define OUTFLOW_OUTPUT                    NO
#define OUTFLOW_OUTPUT_RATE               1.0
#define CLOUD_OUTPUT                      NO
#define CLOUD_OUTPUT_RATE                 1.0


/* --- Not usually changed ---- */
#define MU_CALC                           MU_ANALYTIC
#define CLOUD_TCRIT                       3.0e4
#define CLOUD_MUCRIT                      0.6212407755077543
#define JD_MODE                           JD_GRAD
#define BH_POT_SMOOTH                     4.0

