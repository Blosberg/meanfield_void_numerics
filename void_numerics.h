/*  void_numerics.h--class def's and func decl's (func definitions are in another .h file)
// ---last updated on  Sun Apr 13 23:34:52 CEST 2014  by  bren  at location  bren-Desktop

//  changes from  Sun Apr 13 23:34:52 CEST 2014 : modified the script for explicitly-sized particles (i.e. suitable for dimers, trimers, etc., no longer just coarse-grained full-sized Nucl's.)
*/


 #ifndef __v_numerics  //---check whether it's been defined already so we don't do it twice. 
 #define __v_numerics 

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <sstream>
#include <sys/dir.h>

// #ifndef void_numerics  //---check whether it's been defined already so we don't do it twice. 
// #define void_numerics

using namespace std;

const bool irreversible=false;


//******************************************************************
int space_tpoints_linear(const double tmin, const double tmax, const int nbins, double * tx);
int space_tpoints_logarithmic(const double t0, const double tf, const int nbins, double * tx);
string convertDouble(double number);
string convertInt(double number);

//int  charlength=400; //--the number of characters in the string for our path.
/*********************************************************************************/
class ODEdat
{
public:
// double  *V; //-----------this will be the array of voids.

double * v2; 		//----- v2[x]       = the energetic potential of two particles with 'x' empty sites between them.
double * Bzman_v2;	//----- Bzman_v2[x] = exp(-v2[x]), the boltzmann factor thereof.
int size_v2;		// the number of points in the v2 potential (beyond which it is assumed v2=0)
double * xcoarse;

double rho;		//---the average density.
double maxrho;		//--- the maximum density observed within a single run.

double rhostar;		//---the normalization factor.
double coverage;	//---between [0,1], the fraction of lattice sites covered by a particle.

//! double rhodot;		//---the derivative
double empty_space;
double C1,C2;

double mean;
double std_dev;

double rm; // = r-'minus' -the off rate.
double rp; // = r-'plus'  -the on rate.
double r;  // = r+/r-
double muN;

int a; //---the range of interaction.
int L; //---the size of the system /max void.

double E0;
double t;
double t0, t1; //--the initial step time, and the final simulation time.

bool* printq;		//----the queue of plots to make.
double* tpoints_vdist;	//----the array of points in time to make this plot.
int     total_obs_vdist;	//----the number of plots of the void distribution we should make.
bool  shouldplotvoiddist; //---whether or not we should plot the void distribution at various times.
bool  shouldplotrhos;
int   plotnum; //----the current index of the void-distribution plot.

char cpath;   // [charlength];

int attempt;
string pathout;

bool HNG, SNG, LNG; //--the type of 2-body interaction.

bool boltzmann_on_uphill ;
bool boltzmann_on_add    ;
bool boltzmann_on_removal; //-----this condition is read in from file.

ofstream *log;
int step;
//-------------  functions  ----------------
ODEdat(const double * rates_times, const  int * sizes, const  bool * B, const string pathout);
~ODEdat();
ODEdat& operator= (const ODEdat& param);

int    initialize_V(const int init, double * V);
bool   printtime(void);
int    printout_voiddist(const double * V);
int    get_mean_stddev(const double * V);
double get_empty_space(const double * V);	//--- the fraction of space unoccupied.
double get_rho_anal(const double * V);		//--- the density per system length
double get_C(const double * V);			//--- allocates the two terms of the conserved quantity
double get_rhostar(const double * V);		//--- the normalization factor on the convolution term to impose conservation.

//! double get_rhodot_anal(const double * V);

// bool   printtime(void);
};
//*********************************************************************************

//int time_go(ODEdat *P, ofstream *foutmain, const int init);
int time_go(ODEdat *P, ofstream *foutmain);

int func_HNG (double t, const double V[], double f[], void *params);
int func_SLNG (double t, const double V[], double f[], void *params);

//-----prints out the void distribution as a function of size.

#endif //--this ends the clause as to whether or not this symbol (i.e. this file) has been defined already.
