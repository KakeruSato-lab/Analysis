/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Provide basic functionality for reading input data files.

  Collects a number of functions for opening, reading and assigning 
  initial conditions from user-supplied data file(s).
  The geometry and dimensions of the input grid can be different from 
  the actual grid employed by PLUTO, as long as the coordinate geometry
  transformation has been implemented.
  The input grid and data files should employ the same format and
  conventions employed by PLUTO. 
  
  - Gridfile: coordinates should be written using the PLUTO 4.0 grid format.
  - Datafile: variables should be written in sequence in a single binary 
              file using single or double precision. 
              The file extension must be ".flt" or ".dbl" for the former and
              the latter, respectively.
     
    Note that not all of the variables should be present and the input
    array ::get_var specifies which ones are to be searched for.

  The InputDataSet() initialize the module and by assigning values to 
  global variables such as size, geometry and dimensions of the input grid.
  Data values are read through the function InputDataRead() while
  InputDataInterpolate() can be finally used to map input data onto the
  grid employed by PLUTO using bi- or tri-linear interpolation to fill the 
  data array at the desired coordinate location.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos 
  \date   Aug 27, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* AYW -- 2012-11-15 15:35 JST
 * Reduce to reduce memory requirements. Need one extra for the -1 
 * after the last variable ID in get_var.
 * 2,   if   only rho
 * 3-5, if   rho, vx1{, vx2, vx3}
 */
//#define ID_MAX_NVAR 256
#define ID_MAX_NVAR (1 + (CLOUD_VELOCITY != NONE ? COMPONENTS : 0) + 1)
/* -- AYW */

/* AYW -- 2012-11-29 12:20 JST
 * To use my geometry macros */
#include "pluto_usr.h"

/* -- AYW */

/* AYW -- 2012-11-29 12:20 JST
 * To store input data extents */
double g_idBoxBeg[3];  /**< Lower limits of the input data domain. */
double g_idBoxEnd[3];  /**< Upper limits of the input data domain. */
double g_idnx1, g_idnx2, g_idnx3;
int g_idnvar, *g_idvarindx; /**< Number of variables */
/* -- AYW */


static int id_nvar; /**< Number of variables to be read on input. */
static int id_var_indx[ID_MAX_NVAR]; /**< The variable index. */
static int id_nx1; /**< Size of input grid in the x1 direction. */
static int id_nx2; /**< Size of input grid in the x2 direction. */
static int id_nx3; /**< Size of input grid in the x3 direction. */

static int id_geom;   /**< Geometry of the input grid. */

static double *id_x1; /**< Array of point coordinates of the x1 input grid. */
static double *id_x2; /**< Array of point coordinates of the x2 input grid. */
static double *id_x3; /**< Array of point coordinates of the x3 input grid. */

static double ***Vin[ID_MAX_NVAR]; /**< An array of 3D data values containing the
                                        initial data file variables. */


/* This function shifts a coordinate by a multiple of the input cube width
   DM 11Aug15: Fixed bug. */
double shift_idx(double x, const double xmin, const double xmax, 
                 const double id_x_beg, const double id_x_end) {

    if (x < id_x_beg || x > id_x_end) {
        if (x > xmin && x < xmax) {

            double delta, deltax, intpart, fractpart;
            delta = id_x_end - id_x_beg;
            deltax = x - (x < id_x_beg ? id_x_beg : id_x_end);
            fractpart = modf(fabs(deltax) / delta, &intpart);
            x = (x < id_x_beg ? x + (intpart + 1.0) * delta : x - (intpart + 1.0) * delta);
        }

        if (x < xmin) x = xmin;
        if (x > xmax) x = xmax;
    }

    return x;

}

//--DM--//


/* ********************************************************************* */
void InputDataSet (char *grid_fname, int *get_var)
/*!
 * Initialize access to input data file by assigning values to 
 * grid-related information (geometry, number of points, etc...).
 * This function should be called just once for input-data initialization.
 *
 * \param [in] gname the grid file name
 * \param [in] get_var an array of integers specifying which variables
 *                     have to be read from the input data. 
 * \return Thi function has no return value.
 *
 * The following tasks are performed.
 *********************************************************************** */
{
    int i, ip, nv, success;
    size_t dsize, dcount;
    char *sub_str, sline[256];
    const char delimiters[] = " \t\r\f\n";
    double xl, xr;
    fpos_t file_pos;
    FILE *fp;

#if CLOUDS_MULTI != YES
    print("> Input data:\n\n");
#endif

/* --------------------------------------------------------------------- */
/*! - Scan grid data file and try to determine the grid geometry 
      (::id_geom). Search for tag "GEOMETRY:" and read the word that
      follows.                                                           */
/* --------------------------------------------------------------------- */

    fp = fopen(grid_fname, "r");
    if (fp == NULL) {
        print("! InputDataSet: grid file %s not found\n", grid_fname);
        QUIT_PLUTO(1);
    }
    success = 0;
    while (!success) {
        fgets(sline, 512, fp);
        sub_str = strtok(sline, delimiters);
        while (sub_str != NULL) {
            if (!strcmp(sub_str, "GEOMETRY:")) {
                sub_str = strtok(NULL, delimiters);
                success = 1;
                break;
            }
            sub_str = strtok(NULL, delimiters);
        }
    }

    if (!strcmp(sub_str, "CARTESIAN")) id_geom = CARTESIAN;
    else if (!strcmp(sub_str, "CYLINDRICAL")) id_geom = CYLINDRICAL;
    else if (!strcmp(sub_str, "POLAR")) id_geom = POLAR;
    else if (!strcmp(sub_str, "SPHERICAL")) id_geom = SPHERICAL;
    else {
        print("! InputDataSet: unknown geometry\n");
        QUIT_PLUTO(1);
    }

#if CLOUDS_MULTI != YES
    print("  Input grid file:       %s\n", grid_fname);
    print("  Input grid geometry:   %s\n", sub_str);
#endif

/* --------------------------------------------------------------------- */
/*! - Move file pointer until the first line that does not
      begin with a "#".                                                  */
/* --------------------------------------------------------------------- */

    success = 0;
    while (!success) {
        fgetpos(fp, &file_pos);
        fgets(sline, 512, fp);
        if (sline[0] != '#') success = 1;
    }

    fsetpos(fp, &file_pos);

/* --------------------------------------------------------------------- */
/*! - Start reading number of points and grid coordinates. For the
      input x1 direction these are stored inside the module variables
      ::id_nx1 and ::id_x1.                                              */
/* --------------------------------------------------------------------- */

    /* AYW -- 2013-04-19 12:08 JST
    * Store read data extents in global arrays */

    fscanf(fp, "%d \n", &id_nx1);
    id_x1 = ARRAY_1D(id_nx1, double);
    for (i = 0; i < id_nx1; i++) {
        fscanf(fp, "%d  %lf %lf\n", &ip, &xl, &xr);
        id_x1[i] = 0.5 * (xl + xr);
        if (i == 0)          g_idBoxBeg[IDIR] = xl;   // AYW added
        if (i == id_nx1 - 1) g_idBoxEnd[IDIR] = xr;   // AYW added
    }

    fscanf(fp, "%d \n", &id_nx2);
    id_x2 = ARRAY_1D(id_nx2, double);
    for (i = 0; i < id_nx2; i++) {
        fscanf(fp, "%d  %lf %lf\n", &ip, &xl, &xr);
        id_x2[i] = 0.5 * (xl + xr);
        if (i == 0)          g_idBoxBeg[JDIR] = xl;   // AYW added
        if (i == id_nx2 - 1) g_idBoxEnd[JDIR] = xr;   // AYW added
    }

    fscanf(fp, "%d \n", &id_nx3);
    id_x3 = ARRAY_1D(id_nx3, double);
    for (i = 0; i < id_nx3; i++) {
        fscanf(fp, "%d  %lf %lf\n", &ip, &xl, &xr);
        id_x3[i] = 0.5 * (xl + xr);
        if (i == 0)          g_idBoxBeg[KDIR] = xl;   // AYW added
        if (i == id_nx3 - 1) g_idBoxEnd[KDIR] = xr;   // AYW added
    }
    fclose(fp);

/* -- reset grid with 1 point -- */

    if (id_nx1 == 1) id_x1[0] = 0.0;
    if (id_nx2 == 1) id_x2[0] = 0.0;
    if (id_nx3 == 1) id_x3[0] = 0.0;

    /* AYW -- 2013-04-19 12:08 JST
    * Store read data number of points in global vars */
    g_idnx1 = id_nx1;
    g_idnx2 = id_nx2;
    g_idnx3 = id_nx3;
    /* -- AYW */

    /* DM -- */
    /* In case of CLOUDS_MULTI, overwrite these with input parameters */
#if CLOUDS_MULTI == YES
    g_idBoxBeg[IDIR] = g_inputParam[PAR_WX1L]; g_idBoxEnd[IDIR] = g_inputParam[PAR_WX1H];
    g_idBoxBeg[JDIR] = g_inputParam[PAR_WX2L]; g_idBoxEnd[JDIR] = g_inputParam[PAR_WX2H];
    g_idBoxBeg[KDIR] = g_inputParam[PAR_WX3L]; g_idBoxEnd[KDIR] = g_inputParam[PAR_WX3H];
    /* -- DM */
#endif


#if CLOUDS_MULTI != YES
    print("  Input grid extension:  x1 = [%12.3e, %12.3e] (%d points)\n",
           id_x1[0], id_x1[id_nx1 - 1], id_nx1);
    print("\t\t\t x2 = [%12.3e, %12.3e] (%d points)\n",
           id_x2[0], id_x2[id_nx2 - 1], id_nx2);
    print("\t\t\t x3 = [%12.3e, %12.3e] (%d points)\n",
           id_x3[0], id_x3[id_nx3 - 1], id_nx3);
#endif

/* --------------------------------------------------------------------- */
/*! - Find out how many and which variables we have to read (:id_nvar 
      and ::id_var_indx). 
      Stop counting variables as soon as the first occurrence of "-1" 
      in get_var is encountered                                          */
/* --------------------------------------------------------------------- */

    id_nvar = 0;
    for (nv = 0; nv < ID_MAX_NVAR; nv++) {
        if (get_var[nv] != -1) {
            id_nvar++;
            id_var_indx[nv] = get_var[nv];
        } else {
            break;
        }
    }

#if CLOUDS_MULTI != YES
    print("  Number of variables:   %d\n", id_nvar);
#endif

    /* AYW -- 2014-06-03 20:08 JST
     * Store number of variables and indices in global variables */
    g_idnvar = id_nvar;
    g_idvarindx = ARRAY_1D(ID_MAX_NVAR, int);
    for (nv = 0; nv < id_nvar; nv++) g_idvarindx[nv] = id_var_indx[nv];
    /* -- AYW */

}

/* ********************************************************************* */
void InputDataRead (char *data_fname, char *endianity)
/*!
 * Read input data file and store the contents into the local storage
 * array ::Vin. Memory allocation is also done here.  
 * The grid size and number of variables must have 
 * previously set by calling InputDataSet().
 * 
 * \param [in] data_fname the data file name
 * \param [in] endianity  an input string ("little" or "big") giving 
 *                        the byte-order of how the input data file 
 *                        was originally written.
 *                        If an empty string is supplied, no change is 
 *                        made.
 * \return This function has no return value.
 *********************************************************************** */
{
    int i, j, k, nv, swap_endian = NO;
    size_t dsize, dcount;
    double udbl;
    float uflt;
    char ext[] = "   ";
    FILE *fp;

/* ----------------------------------------------------
             Check endianity 
   ---------------------------------------------------- */

    if ((!strcmp(endianity, "big") && IsLittleEndian()) ||
        (!strcmp(endianity, "little") && !IsLittleEndian())) {
        swap_endian = YES;
    }

#if CLOUDS_MULTI != YES
    print("  Input data file:       %s (endianity: %s) \n",
           data_fname, endianity);
#endif

/* ------------------------------------------------------
    Get data type from file extensions (dbl or flt).
   ------------------------------------------------------ */

    dcount = strlen(data_fname);
    for (i = 0; i < 3; i++) ext[i] = data_fname[dcount - 3 + i];

    if (!strcmp(ext, "dbl")) {
#if CLOUDS_MULTI != YES
        print("  Precision:             (double)\n");
#endif
        dsize = sizeof(double);
    } else if (!strcmp(ext, "flt")) {
#if CLOUDS_MULTI != YES
        print("  Precision:             (single)\n");
#endif
        dsize = sizeof(float);
    } else {
        print("! InputDataRead: unsupported data type '%s'\n", ext);
        QUIT_PLUTO(1);
    }

/* -------------------------------------------------------
     Read and store data values
   ------------------------------------------------------- */

    fp = fopen(data_fname, "rb");
    if (fp == NULL) {
        print("! InputDataRead: file %s does not exist\n", data_fname);
        QUIT_PLUTO(1);
    }
    for (nv = 0; nv < id_nvar; nv++) {
        if (Vin[nv] == NULL) Vin[nv] = ARRAY_3D(id_nx3, id_nx2, id_nx1, double);

        dcount = 1;


        if (dsize == sizeof(double)) {
            for (k = 0; k < id_nx3; k++) {
                for (j = 0; j < id_nx2; j++) {
                    for (i = 0; i < id_nx1; i++) {
                        if (fread(&udbl, dsize, dcount, fp) != dcount) {
                            print("! InputDataRead: error reading data %d.\n", nv);
                            break;
                        }
                        if (swap_endian) SWAP_VAR(udbl);
                        Vin[nv][k][j][i] = udbl;
                    }
                }
            }
        } else {
            for (k = 0; k < id_nx3; k++) {
                for (j = 0; j < id_nx2; j++) {
                    for (i = 0; i < id_nx1; i++) {
                        if (fread(&uflt, dsize, dcount, fp) != dcount) {
                            print("! InputDataRead: error reading data %d.\n", nv);
                            break;
                        }
                        if (swap_endian) SWAP_VAR(uflt);
                        Vin[nv][k][j][i] = uflt;
                    }
                }
            }
        }
    }

    fclose(fp);

#if CLOUDS_MULTI != YES
    print("\n");
#endif
}

/* ********************************************************************* */
void InputDataInterpolate (double *vs, double x1, double x2, double x3)
/*!
 * Perform bi- or tri-linear interpolation on external
 * dataset to compute vs[] at the given point {x1,x2,x3}.
 *
 * \param [in] vs interpolated value
 * \param [in] x1 coordinate point at which at interpolates are desired
 * \param [in] x2 coordinate point at which at interpolates are desired
 * \param [in] x3 coordinate point at which at interpolates are desired
 * \return This function has no return value.
 *
 * The function performs the following tasks. 
 *********************************************************************** */
{
    int il = 0, jl = 0, kl = 0;
    int ih, jh, kh;
    int im, jm, km;
    int i, j, k, nv, inv;
    double xx, yy, zz;
    double ***V;
    /* AYW -- 2012-11-29 11:53 JST
     * Do geometric velocity conversions here */
    double vel, lorentz;
    int do_vel = 0;
    /* -- AYW */

    double deltax2, deltax3, deltax1, fractpart, intpart;



// TODO: Actually do this check
/* --------------------------------------------------------------------- */
/*! - Convert PLUTO coordinates to input grid geometry if necessary.     */
/* AYW --  2012-11-15 15:48 JST
 * Could use my generic macros for coordinates here.
 * Good check as to whether they're correct.
 * -- AYW */
/* --------------------------------------------------------------------- */


/* AYW -- 2012-11-29 12:33 JST
 * My version. */
    if (id_geom == CARTESIAN) {
        double xc1, xc2, xc3;
        xc1 = CART1(x1, x2, x3);
        xc2 = CART2(x1, x2, x3);
        xc3 = CART3(x1, x2, x3);
        x1 = xc1;
        x2 = xc2;
        x3 = xc3;
    }

    else if (id_geom == CYLINDRICAL) {
        double R, z, phi;
        R = CYL1(x1, x2, x3);
        phi = CYL2(x1, x2, x3);
        z = x3;
        x1 = R;
        x2 = z;
        x3 = phi;
    }

    else if (id_geom == POLAR) {
        double R, phi, z;
        R = POL1(x1, x2, x3);
        phi = POL2(x1, x2, x3);
        z = POL3(x1, x2, x3);
        x1 = R;
        x2 = phi;
        x3 = z;
    }

    else if (id_geom == SPHERICAL) {
        double R, theta, phi;
        R = SPH1(x1, x2, x3);
        theta = SPH2(x1, x2, x3);
        phi = SPH3(x1, x2, x3);
        x1 = R;
        x2 = theta;
        x3 = phi;
    }

/* -- AYW */

//
//  #if GEOMETRY == CARTESIAN
//   if (id_geom == GEOMETRY) {
//
//     /* same coordinate system: nothing to do */
//
//   }else if (id_geom == CYLINDRICAL) {
//     double R, z, phi;
//     R   = sqrt(x1*x1 + x2*x2);
//     phi = atan2(x2,x1);
//     if (phi < 0.0) phi += 2.0*CONST_PI;
//     z   = x3;
//
//     x1 = R; x2 = z; x3 = phi;
//   }else if (id_geom == POLAR) {
//     double R, phi, z;
//     R   = sqrt(x1*x1 + x2*x2);
//     phi = atan2(x2,x1);
//     if (phi < 0.0) phi += 2.0*CONST_PI;
//     z   = x3;
//
//     x1 = R; x2 = phi; x3 = z;
//   }else if (id_geom == SPHERICAL){
//     double r, theta, phi;
//     r     = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
//     r     = sqrt(r);
//     theta = acos(x3/r);
//     phi   = atan2(x2,x1);
//     if (phi   < 0.0) phi   += 2.0*CONST_PI;
//     if (theta < 0.0) theta += 2.0*CONST_PI;
//
//     x1 = r; x2 = theta; x3 = phi;
//   }else{
//     print ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
//     QUIT_PLUTO(1);
//   }
//  #elif GEOMETRY == CYLINDRICAL
//   if (id_geom == GEOMETRY) {
//
//     /* same coordinate system: nothing to do */
//
//   }else if (id_geom == SPHERICAL) {
//     double r, theta, phi;
//     r     = D_EXPAND(x1*x1, + x2*x2, + 0.0);
//     r     = sqrt(r);
//     theta = acos(x2/r);
//     phi   = 0.0;
//     if (theta < 0.0) theta += 2.0*CONST_PI;
//
//     x1 = r; x2 = theta; x3 = phi;
//   }else{
//     print ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
//     QUIT_PLUTO(1);
//   }
//  #elif GEOMETRY == POLAR
//   if (id_geom == GEOMETRY) {
//
//     /* same coordinate system: nothing to do */
//
//   }else if (id_geom == CARTESIAN) {
//     double x, y, z;
//     x = x1*cos(x2);
//     y = x1*sin(x2);
//     z = x3;
//
//     x1 = x; x2 = y; x3 = z;
//   }else{
//     print ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
//     QUIT_PLUTO(1);
//   }
//  #elif GEOMETRY == SPHERICAL
//   if (id_geom == GEOMETRY) {
//
//     /* same coordinate system: nothing to do */
//
//
//   }else if (id_geom == CARTESIAN) {
//     double x, y, z;
//     x = x1*sin(x2)*cos(x3);
//     y = x1*sin(x2)*sin(x3);
//     z = x1*cos(x2);
//
//     x1 = x; x2 = y; x3 = z;
//   }else{
//     print ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
//     QUIT_PLUTO(1);
//   }
//  #endif

/* --------------------------------------------------------------------- */
/*! - Make sure point (x1,x2,x3) does not fall outside input grid range. 
      Limit to input grid edge otherwise.                                */
/* --------------------------------------------------------------------- */


#if (CLOUD_REPEAT == YES)
    /*
     DM 23 feb, 2015: Check if PLUTO grid point is outside input grid.
    If yes, then check if PLUTO grid is outside maximum zone of clouds
    defined in definitions_usr.h. If no, then fold PLUTO grid back on
    to appropriate location in input grid, assuming periodic BC.

    XL, XR are domain limits of input grid

    if x < XL
    |-----------|------------|--------------|
          x     XL           XR
          -------------->
    x' = x + (int((x - xL) / Delta) + 1) * Delta , Delta = XR - XL

    if x > XR
    |-----------|------------|--------------|
                XL           XR        x
                        <--------------
    x' = x - (int((x - xR) / Delta) + 1) * Delta

    Rationale: Translate input grid by n steps such that PLUTO grid point
    falls within translated input grid limits
    i.e. XR' = XR + n * Delta, n = int((x - XR) / Delta) + 1
    New x' = XR - (XR' - x) = x - n * Delta

    Same for lower boundary, results in a difference in sign.

    *** To fix: if PLUTO grid is larger than max domain defined in definitions_usr.h
    then PLUTO grid is given value eq to boundary value of input grid.
    */

    /* DM -- */
      D_EXPAND(x1 = shift_idx(x1, g_inputParam[PAR_WX1L], g_inputParam[PAR_WX1H], id_x1[0], id_x1[id_nx1 - 1]);,
               x2 = shift_idx(x2, g_inputParam[PAR_WX2L], g_inputParam[PAR_WX2H], id_x2[0], id_x2[id_nx2 - 1]);,
               x3 = shift_idx(x3, g_inputParam[PAR_WX3L], g_inputParam[PAR_WX3H], id_x3[0], id_x3[id_nx3 - 1]););
    /* -- DM */

#else
    /* Default PLUTO */
    D_EXPAND(if (x1 < id_x1[0]) x1 = id_x1[0];
             else if (x1 > id_x1[id_nx1 - 1]) x1 = id_x1[id_nx1 - 1];,

             if (x2 < id_x2[0]) x2 = id_x2[0];
             else if (x2 > id_x2[id_nx2 - 1]) x2 = id_x2[id_nx2 - 1];,

             if (x3 < id_x3[0]) x3 = id_x3[0];
             else if (x3 > id_x3[id_nx3 - 1]) x3 = id_x3[id_nx3 - 1];)
#endif


/* --------------------------------------------------------------------- */
/*! - Use table lookup by binary search to  find the indices 
      il, jl and kl such that grid points of PLUTO fall between 
      [il, il+1], [jl, jl+1], [kl, kl+1].                                */
/* --------------------------------------------------------------------- */

    il = 0;
    ih = id_nx1 - 1;
    while (il != (ih - 1)) {
        im = (il + ih) / 2;
        if (x1 <= id_x1[im]) ih = im;
        else il = im;
    }

    if (id_nx2 > 1) {
        jl = 0;
        jh = id_nx2 - 1;
        while (jl != (jh - 1)) {
            jm = (jl + jh) / 2;
            if (x2 <= id_x2[jm]) jh = jm;
            else jl = jm;
        }
    }

    if (id_nx3 > 1) {
        kl = 0;
        kh = id_nx3 - 1;
        while (kl != (kh - 1)) {
            km = (kl + kh) / 2;
            if (x3 <= id_x3[km]) kh = km;
            else kl = km;
        }
    }

/* --------------------------------------------------------------------- */
/*! - Define normalized coordinates between [0,1]:
      - x[il+1] < x1[i] < x[il+1] ==> 0 < xx < 1
      - y[jl+1] < x2[j] < y[jl+1] ==> 0 < yy < 1
      - z[kl+1] < x3[k] < z[kl+1] ==> 0 < zz < 1                            */
/* --------------------------------------------------------------------- */

    xx = yy = zz = 0.0; /* initialize normalized coordinates */

    if (id_nx1 > 1) xx = (x1 - id_x1[il]) / (id_x1[il + 1] - id_x1[il]);
    if (id_nx2 > 1) yy = (x2 - id_x2[jl]) / (id_x2[jl + 1] - id_x2[jl]);
    if (id_nx3 > 1) zz = (x3 - id_x3[kl]) / (id_x3[kl + 1] - id_x3[kl]);

/* --------------------------------------------------------------------- */
/*! - Perform bi- or tri-linear interpolation.                           */
/* --------------------------------------------------------------------- */

    for (nv = 0; nv < id_nvar; nv++) {
        inv = id_var_indx[nv];

        /* AYW -- 2012-11-29 11:53 JST
         * Check if velocity geometric conversion needs to be done. */
        if ((inv >= VX1) && (inv <= VX3)) do_vel = 1;
        /* --AYW */

        V = Vin[nv];
        vs[inv] = V[kl][jl][il] * (1.0 - xx) * (1.0 - yy) * (1.0 - zz)
                  + V[kl][jl][il + 1] * xx * (1.0 - yy) * (1.0 - zz);
        if (id_nx2 > 1) {
            vs[inv] += V[kl][jl + 1][il] * (1.0 - xx) * yy * (1.0 - zz)
                       + V[kl][jl + 1][il + 1] * xx * yy * (1.0 - zz);
        }
        if (id_nx3 > 1) {
            vs[inv] += V[kl + 1][jl][il] * (1.0 - xx) * (1.0 - yy) * zz
                       + V[kl + 1][jl][il + 1] * xx * (1.0 - yy) * zz
                       + V[kl + 1][jl + 1][il] * (1.0 - xx) * yy * zz
                       + V[kl + 1][jl + 1][il + 1] * xx * yy * zz;
        }
    }

/* AYW -- 2012-11-29 11:53 JST 
 * Do geometric velocity conversions here. Convert from
 * data geometry back to PLUTO geometry 
 * NOTE: Assumptions on alignment made.
*/
    if (do_vel) {

        if (id_geom == CARTESIAN) {
            double v1, v2, v3;
            v1 = vs[VX1];
            v2 = vs[VX2];
            v3 = vs[VX3];
            EXPAND(vs[VX1] = VCART_1(x1, x2, x3, v1, v2, v3);,
                   vs[VX2] = VCART_2(x1, x2, x3, v1, v2, v3);,
                   vs[VX3] = VCART_3(x1, x2, x3, v1, v2, v3););
        }

        else if (id_geom == CYLINDRICAL) {
            double v1, v2, v3;
            v1 = vs[VX1];
            v2 = vs[VX2];
            v3 = vs[VX3];
            EXPAND(vs[VX1] = VCYL_1(x1, x2, x3, v1, v2, v3);,
                   vs[VX2] = VCYL_2(x1, x2, x3, v1, v2, v3);,
                   vs[VX3] = VCYL_3(x1, x2, x3, v1, v2, v3););
        }

        else if (id_geom == POLAR) {
            double v1, v2, v3;
            v1 = vs[VX1];
            v2 = vs[VX2];
            v3 = vs[VX3];
            EXPAND(vs[VX1] = VPOL_1(x1, x2, x3, v1, v2, v3);,
                   vs[VX2] = VPOL_2(x1, x2, x3, v1, v2, v3);,
                   vs[VX3] = VPOL_3(x1, x2, x3, v1, v2, v3););
        }

        else if (id_geom == SPHERICAL) {
            double v1, v2, v3;
            v1 = vs[VX1];
            v2 = vs[VX2];
            v3 = vs[VX3];
            EXPAND(vs[VX1] = VSPH_1(x1, x2, x3, v1, v2, v3);,
                   vs[VX2] = VSPH_2(x1, x2, x3, v1, v2, v3);,
                   vs[VX3] = VSPH_3(x1, x2, x3, v1, v2, v3););
        }

        else {
            print("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
            QUIT_PLUTO(1);
        }


    }
    /* -- AYW */


}

/* ********************************************************************* */
void InputDataFree (void)
/*!
 * Free memory stored by user-supplied data.
 *
 *********************************************************************** */
{
    int nv;
    for (nv = 0; nv < id_nvar; nv++) {
        free((char *) Vin[nv][0][0]);
        free((char *) Vin[nv][0]);
        free((char *) Vin[nv]);
    }
}

