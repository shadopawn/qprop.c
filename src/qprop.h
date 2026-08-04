/*******************************************************************************
    qprop.c: a simple and lightweight library for propeller aerodynamic analysis

    The blade-element/momentum equations are solved with the method of
    S. Andrew Ning (Wind Energy 17(9), 2014): a bracketed one-dimensional root
    find on the local inflow angle, with momentum theory switched between the
    momentum, empirical (Buhl) and reversed-disk-flow branches. Well-suited for
    rotors that operate at low Reynolds numbers and do not feature complex 3D
    effects; handles propeller, turbine/windmill and descent windmill brake
    operation.

    Key characteristics:
    - Lightweight and portable: contained in a single file with no dependencies
    - No file I/O required: perform analyses and retrieve results without
      writing input files or reading output files
    - Easy to use: simply copy the library into your project directory
    - Accurate airfoil polars: unlike the original QPROP, which requires users
      to tune oversimplified analytic models, qprop.c gets the aerodynamic
      coefficients of the airfoils by interpolating XFoil polars

    This file is a header for the precompiled shared libraries,
    replicating the definitions in qprop.c.

    Author: Andrea Pavan
    License: MIT
*******************************************************************************/


//---------------------
//  DATA STRUCTURES
//---------------------

//data structure for polars
typedef struct {
    double Re;          //Reynolds number
    double* alpha;      //array of angle of attacks (rad)
    double* CL;         //array of lift coefficients - same size as alpha
    double* CD;         //array of drag coefficients - same size as alpha
    int size;           //number of points in the polar
    //Uniform-alpha resampling of the polar, built once at load time. Lookups
    //then cost index arithmetic instead of a linear scan over alpha. The grid
    //is at least as fine as the source data, so it reproduces the scanned
    //interpolation to within rounding inside the data range, and it also bakes
    //in the post-stall extrapolation outside it. These fields are appended so
    //that the offsets of the ones above are unchanged.
    double* gridCL;     //NULL when the grid has not been built
    double* gridCD;
    int gridN;          //number of grid points
    double gridA0;      //alpha of the first grid point (rad)
    double gridDA;      //grid spacing (rad)
} Polar;

//data structure for airfoils
typedef struct {
    Polar** polars;     //array of pointers to polars - typically at different Re
    int size;           //number of polars in the airfoil
} Airfoil;

//data structure for blade sections
typedef struct {
    double c;           //chord length (m)
    double beta;        //twist angle (rad)
    double r;           //radial distance (m)
    Airfoil airfoil;   //local airfoil data
} Section;

//data structure for warm-start state
//
//Holds the converged inflow angle of every blade element from the previous
//solve, so that a following solve at a nearby operating point can seed a secant
//iteration from it instead of bracketing from scratch. Intended for real-time
//use, where the operating point moves very little between frames: measured at
//1000 Hz on a quadcopter trajectory it cuts the solve from about 23 residual
//evaluations per element to about 4, a little over five times faster.
//
//A rotor has no state until prop_state_new() is called on it, so warm starting
//is opt-in and the default solve stays stateless and reproducible.
//
//With a state attached the solver becomes path dependent: near a bifurcation,
//where both halves of the phi domain admit a root, the answer depends on how
//the operating point was reached. That is not unphysical - real rotors show
//hysteresis entering and leaving the vortex ring state - but the result is no
//longer a function of the operating point alone, so anything that needs to be
//reproducible (a replay, a regression test) must checkpoint the state too.
typedef struct {
    double* phi;        //converged inflow angle per element (rad)
    int nelems;         //number of elements this state was sized for
    int valid;          //0 until a first solve has filled phi
} PropState;

//data structure for rotors
typedef struct {
    double D;           //rotor diameter (m)
    int B;              //number of blades
    int nsections;      //number of sections discretizing a blade
    Section* sections;  //array of sections discretizing a blade
    //Warm-start state, owned by the rotor and freed with it. NULL unless
    //prop_state_new() has been called, in which case qprop() seeds each
    //element from the previous solve. Appended so that the offsets of the
    //fields above are unchanged.
    PropState* state;
} Rotor;

//data structure for qprop output
typedef struct {
    double T;           //overall thrust (N)
    double Q;           //overall torque (N-m)
    double CT;          //thrust coefficient
    double CP;          //power coefficient
    double J;           //advance ratio
    double* residuals;  //array of elements residuals
    double* Gamma;      //array of elements circulations
    double* lambdaw;    //array of local wake advance ratios
    double* r;          //array of elements radial distance (m)
    double* W;          //array of local velocities (m/s)
    double* phi;        //array of local inflow angle (rad)
    double* dTdr;       //array for blade thrust distribution (N/m)
    double* dQdr;       //array for blade torque distribution (N-m/m)
    int nelems;         //number of elements discretizing a blade
} RotorPerformance;


//-----------------
//  FREE MEMORY
//-----------------

//FREE_POLAR frees the memory allocated in a Polar structure
//Input:
//  - currentpolar (Polar*): pointer to a polar that is no longer needed
//Output:
//  - none
void free_polar(Polar* currentpolar);

//FREE_AIRFOIL frees the memory allocated in an Airfoil structure
//Input:
//  - currentairfoil (Airfoil*): pointer to an airfoil that is no longer needed
//Output:
//  - none
void free_airfoil(Airfoil* currentairfoil);

//FREE_ROTOR frees the memory allocated in a Rotor structure
//Input:
//  - currentrotor (Rotor*): pointer to a rotor that is no longer needed
//Output:
//  - none
void free_rotor(Rotor* currentrotor);

//FREE_ROTOR_PERFORMANCE frees the memory allocated in a qprop output
//Input:
//  - perf (RotorPerformance*): pointer to a qprop output that is no longer needed
//Output:
//  - none
void free_rotor_performance(RotorPerformance* perf);


//---------------------------
//  FUNCTION DECLARATIONS
//---------------------------

//DEG2RAD converts degrees to radians
//Input:
//  - deg (double): angle in degrees
//Output:
//  - (double): angle in radians
double deg2rad(double deg);

//READ_XFOIL_POLAR_FROM_FILE reads an airfoil polar from a text file
//Input:
//  - filename (array of char): name of the txt file containing the polar data
//Output:
//  - (Polar*): pointer to the polar data
//Notes:
//  - The file is assumed to be in the XFoil/XFLR5 format:
//      - Reynolds number on a line containing "Re =", ignoring spaces
//      - A table of alpha, CL and CD values, ordered by alpha (from min to max)
//      - Alpha values in the first column, CL in the second and CD in the third
//      - No empty lines between values in the table
//  - The file content is not thoroughly checked for errors
//  - This function internally allocates memory for the Polar structure arrays
//    alpha, CL and CD using malloc and realloc.
//    It is the caller's responsibility to free this memory by calling
//    unload_polar_from_memory(Polar*) when it is no longer needed
Polar* read_xfoil_polar_from_file(const char *filename);

//IMPORT_XFOIL_POLARS imports airfoil polars from multiple text files
//Input:
//  - filenames (array of (array of char)): list of files containing polar data
//  - number_of_files (int): number of files in the array
//Output:
//  - (Airfoil*): pointer to the imported airfoil polars
//Notes:
//  - All files are assumed to be in the XFoil/XFLR5 format
//    (see the notes above "read_xfoil_polar_from_file")
//  - Safety checks on user input are not implemented yet
//  - The content of each file is not checked
//  - This function internally allocates memory for the Airfoil structure arrays
//    and for each Polar using malloc and realloc.
//    It is the caller's responsibility to free this memory when it is no longer
//    needed, by calling unload_airfoil_from_memory(Airfoil*)
Airfoil* import_xfoil_polars(const char *filenames[], int number_of_files);

//ANALYTIC_POLAR_CURVES generates polars using the simple analytic model
//described by Drela in the QPROP user guide
//Input:
//  - CL0 (double): zero-lift lift coefficient
//  - CL_a (double): lift curve slope
//  - CLmin (double): minimum lift coefficient
//  - CLmax (double): maximum lift coefficient
//  - CD0 (double): zero-lift drag coefficient
//  - CD2u (double): quadratic coefficient in the drag formula
//  - CD2l (double): quadratic coefficient in the drag formula
//  - CLCD0 (double): lift coefficient at minimum drag
//  - REref (double): reference Reynolds number for all the coefficients above
//  - REexp (double): Reynolds number exponent (suggested: -0.5)
//Output:
//  - (Airfoil*): pointer to generated polar curves
Airfoil* analytic_polar_curves(double CL0, double CL_a, double CLmin, double CLmax,
                               double CD0, double CD2u, double CD2l, double CLCD0,
                               double REref, double REexp);

//IMPORT_ROTOR_GEOMETRY_APC reads a propeller geometry from an APC PE0 file
//Input:
//  - filename (array of char): name of the PE0 file containing the geom data
//  - airfoil (Airfoil*): pointer to an airfoil
//Output:
//  - (Rotor*): pointer to imported rotor geometry with the given airfoil
//Notes:
//  - the file is assumed to be downloaded from the official APC website
Rotor* import_rotor_geometry_apc(const char *filename, Airfoil *airfoil);

//IMPORT_ROTOR_GEOMETRY_UIUC reads a propeller geometry from an UIUC txt file
//Input:
//  - filename (array of char): name of the txt file containing the geom data
//  - airfoil (Airfoil*): pointer to an airfoil
//  - D (double): rotor diameter in m
//  - B (int): number of blades
//Output:
//  - (Rotor*): pointer to imported rotor geometry with the given properties
//Notes:
//  - the file is assumed to be downloaded from the official APC website
Rotor* import_rotor_geometry_uiuc(const char *filename, Airfoil* airfoil, double D, int B);

//REFINE_ROTOR_SECTIONS creates a propeller geometry with the specified number
//of equally-spaced sections
//Input:
//  - oldrotor (Rotor*): pointer to the reference rotor geometry
//  - nsections (int): desired number of sections
//Output:
//  - newrotor (Rotor*): pointer to the new rotor geometry
Rotor* refine_rotor_sections(Rotor* oldrotor, int nsections);

//QPROP runs the QProp algorithm as described by Drela for each blade element
//Input:
//  - rotor (Rotor*): pointer to a rotor
//  - Uinf (double): freestream velocity in m/s
//  - Omega (double): rotor speed in rad/s
//  - tol (double): stopping criterion tolerance (suggested value: 1e-6)
//  - itmax (int): maximum number of iterations (suggested value: 100)
//  - rho (double): air density in kg/m3 (suggested value: 1.225)
//  - mu (double): air dynamic viscosity in Pa-s (suggested value: 1.81e-5)
//  - a (double): speed of sound in m/s (suggested value: 340.0) - set to 0 to disable Mach correction
//Output:
//  - (RotorPerformance*): pointer to the QProp outputs
//Notes:
//  - the current implementation assumes that there is no externally-induced
//    tangential velocity (Ut = 0)
//  - when rotor->state is set (see prop_state_new) each element is seeded from
//    the previous solve and the state is updated in place on the rotor;
//    otherwise the solve is stateless
RotorPerformance* qprop(Rotor* rotor, double Uinf, double Omega, double tol, int itmax, double rho, double mu, double a);

//PROP_STATE_NEW attaches a warm-start state to a rotor, enabling warm starting
//Input:
//  - rotor (Rotor*): the rotor to attach the state to
//Output:
//  - (PropState*): pointer to the attached state, also reachable as
//    rotor->state, or NULL on failure
//Notes:
//  - subsequent calls to qprop() on this rotor seed each element from the
//    previous solve; see the note on PropState about path dependence
//  - the rotor owns the state and free_rotor() releases it, so there is no
//    need to free it separately; calling this twice is harmless
//  - one state per rotor instance: give each propeller its own Rotor
PropState* prop_state_new(Rotor* rotor);

//PROP_STATE_FREE releases a rotor's warm-start state and returns it to
//stateless solving
//Input:
//  - rotor (Rotor*): the rotor whose state should be released
//Output:
//  - none
void prop_state_free(Rotor* rotor);

//PROP_STATE_RESET discards the stored inflow angles without freeing the state,
//so the next solve starts cold and later ones warm again
//Input:
//  - rotor (Rotor*): the rotor whose state should be invalidated
//Output:
//  - none
//Notes:
//  - use after a discontinuity (a propeller stopping, a teleport) where the
//    previous angles are no longer a useful guess
void prop_state_reset(Rotor* rotor);

