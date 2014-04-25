/*--- VOID NUMERICS DRIVER SCRIPT -TESTS THE EQUATIONS DERIVED FOR SNG IN ANALOGY WITH REDNER -----
// ---last updated on  Fri Apr 25 17:08:48 CEST 2014  by  ga79moz  at location  TUM , murter

//  changes from  Fri Apr 25 17:08:48 CEST 2014 : optimized the process by iteratively cutting off large voids that have gone negative. cutoff is in func_SLNG"

//  changes from  Sun Apr 13 23:34:52 CEST 2014 : modified the script for explicitly-sized particles (i.e. suitable for dimers, trimers, etc., no longer just coarse-grained full-sized Nucl's.)

//  changes from  Tue Nov 26 19:18:48 CET 2013 : ran valgrind to clean up mem. management. Added a destructor to the ODEdata class.
//-----------------------------------------------------------------------------------------------*/


#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <dirent.h>



//------------------------------------------
#include "void_numerics.h"
#include "void.h"
#include <bren_lib.h>
//------------------------------------------


using namespace std;

double const pi   = 3.14159265358979323846264;


//******************************************************************

double interpolate_mu_from_rhoi(const double irho_target,const int w, const double E0, const int kHNG, const string NGtype, const string CGF_feedin_string);
//******************************************************************


int main(int argc, char *argv[])
{
int     L=0, a=0, j=0, test_result=0; 
double  t_export=0.0;
//---the size of the system, # of plots to make, size of particle, size of footprint, counting index, dummy.
double  t0, t1, rm, E0,  muN_input, muN;
bool    shouldplotvoiddist, shouldplotrhos, HNG, SNG, LNG, should_check_neg;
bool    boltzmann_on_uphill, boltzmann_on_add, boltzmann_on_removal, should_import_IC, should_export_IC;

string  pathout;

//--'command line parameters' ---------------------------------------------------
if (argc != 4)
	{
	cout << "\n error: expecting 3 command-line argument (1) muN (2) E0, and (3) NG type. Instead we got " << argc << " arguments which were:\n";
	for (j=1;j<argc;j++)
		{
		cout << argv[j] << ",  ";
		}
	cout << "\n exiting\n\n";
	exit(1);
	}
muN_input     = atof(argv[1]); 
E0            = atof(argv[2]); //---atof is for floats and doubles.
string NGtype = argv[3];

//----------------------retrieve input data-------------------------------

int parity_check=-1, charlength=400; //--the number of characters in the string for our path.
char cpath[charlength];
clear_charray(cpath, charlength );

sprintf(cpath, "void_numerics.in");
ifstream fin(cpath);
if(fin.fail())
	{
	cout << "\n ERROR: cannot locate input file.\n";
	exit(1); 
	}

string BZflag; //-----flag to determine which reactions are weighted by boltzmann factors.
//--------------------------NOW READ IN -----------------------------------

fin  >> L;			// size of the system after coarse-graining (therefore size of max void.)
fin  >> t0 >> t1;		// initial time step, termination time.
fin  >> a;			// size of the footprint

//--NB: removed the read-in of E0 and muN_input, they are now command-line inputs.

fin  >> shouldplotvoiddist >> shouldplotrhos >> should_check_neg;		// boolean should we print the void dist profiles.

fin  >> should_import_IC   >> should_export_IC >> t_export;

fin  >> BZflag;
fin  >> parity_check;
fin  >> pathout;
fin.close();

if (parity_check != 88855888)
	{
	cout << "\n ERROR: in ordering of input parameters. exiting. \n"; 
	exit(1);
	}
/*
if ( int(should_export_IC) + int(should_import_IC) !=1)
	{
	cout << "\n ERROR: precisely one of import/export should be true. \n"; 
	exit(1);
	}
*/
//----------------------------------------------------------
string BZ;
if (BZflag =="boltzmann_on_uphill")
	{
	boltzmann_on_uphill  = true;
	boltzmann_on_add     = false;
	boltzmann_on_removal = false;
	BZ="uphill";
	}
else if(BZflag == "boltzmann_on_add")
	{
	boltzmann_on_uphill  = false;
	boltzmann_on_add     = true;
	boltzmann_on_removal = false;
	BZ="add";
	}
else if(BZflag == "boltzmann_on_removal")
	{
	boltzmann_on_uphill  = false;
	boltzmann_on_add     = false;
	boltzmann_on_removal = true;
	BZ="remove";
	}
else
	{
	cout << "\n ERROR: bolzmann criteria is : " << BZflag << " -which is either conflicting or undefined.\n";
	exit(1);
	}

//---------------------------------------------------------

if(irreversible)
	{
	rm=0.0;		// rate off ALWAYS =0
	}
else
	{
	rm=1.0;		// rate off ALWAYS =1 
	}

//---------------------------------------------------------


//---   don't use this anymore: muN is now a direct command-line input. 
//---   muN = muN_input + 1.0*(TASKID-1); 
muN = muN_input; 
//---   this used to be how I would stagger the muN values according to input for parallelization.

//---this is done for both SNG and HNG cases.

int init=L;

if(NGtype=="SNG")
	{
	SNG=true; HNG=false; LNG=false;
	}
else if(NGtype=="LNG")
	{
	LNG=true; SNG=false; HNG=false;
	}
else if(NGtype=="HNG")
	{
	 HNG=true; SNG=false; LNG=false;
	if ( fmod (a,1) > 0.0000001) //----make sure it's an integer
		{
		cout << "\n ERROR: a = " << a << " value is not an integer despite HNG selection. Exiting.\n\n";
		exit(1);
		}
	//---else just continue as normal.
	}
else
	{
	cout << "\n ERROR: ambiguous NG type. exiting. \n"; 
	exit(1);
	}

//---------------------CHECK IF SYSTEM SIZE IS LARGE ENOUGH

if( L < 5*a  )
	{
	cout << "\n WARNING: system size is less than 5 particles. exiting.\n";
//!	cout << " the interaction range of the coarse-grained particles. too small. Exiting.\n";
//!	exit(1);
	}


//-------------------------- setup the log file--------------------------

if ( DirectoryExists( pathout.c_str()  ) )
	{
	cout << "\n output path already exists. \n";
	}
else
	{
	cout << "\n output path does not exist -creating it. \n";
	clear_charray(cpath, charlength );

	sprintf(cpath, "mkdir \"./%s\"",pathout.c_str());
	system(cpath);
	}

clear_charray(cpath, charlength );

sprintf(cpath, "%svoid_job.log",pathout.c_str());

ofstream *flog = new ofstream(cpath);
ofstream v2pot_output;

//----------- DOCUMENTATION ------------------------------------------------------------

*flog << "\n for this run, muN=" << muN;
*flog << "\n muN_input = " << muN_input;
*flog << "\n a = " << a << ", E0 = " << E0;
*flog << "\n before any IC imports, input file specifies L=" << L;
*flog << "\n NGtype = " << NGtype;
*flog << "\n t0 = " << t0 ;
*flog << "\n t1 = " << t1 ;
*flog << "\n shouldplotvoiddist = " << shouldplotvoiddist << ", shouldplotrhos = " << shouldplotrhos;
*flog << "\n BZflag =" << BZflag << endl;

//--------------------- SET UP PARAMETERS --------------------------------------------------

double rates_times[6];
rates_times[0] = t0;
rates_times[1] = t1;
rates_times[2] = rm;
rates_times[3] = muN;
rates_times[4] = E0;
rates_times[5] = t_export;

int sizes[2];
sizes[0] = L;
sizes[1] = a;	//--- the finite size of the particles.
		//--- phy_bound should start off being the same as "L"

bool B[11];
B[0]  = shouldplotvoiddist;
B[1]  = shouldplotrhos;
B[2]  = HNG;
B[3]  = SNG;
B[4]  = LNG;
B[5]  = boltzmann_on_add;
B[6]  = boltzmann_on_removal;
B[7]  = boltzmann_on_uphill;
B[8]  = should_check_neg;
B[9]  = should_import_IC;
B[10] = should_export_IC;


//------------------------------------------------------------------------

ODEdat* P;   //----'P' is simply our catch-all data structure.
ofstream *foutmain;
ofstream *maxrhovals;

if( ! should_export_IC)
	{
	clear_charray(cpath, charlength );
	sprintf(cpath, "%sovershootdat-%sBZ%s_muN_E0_maxrho_finalrho.txt", pathout.c_str(), NGtype.c_str(), BZ.c_str());
	maxrhovals = new ofstream(cpath);	
	}

// ----------------- @@@ HERE IS WHERE WE SHOULD INSERT REPETITION OVER EPSILON SCAN. -----------

P=new ODEdat(rates_times, sizes, B, pathout);   // ---THE 2-BODY INTERACTION IS CALCULATED 
						// ---AND ASSIGNED IN THE CONSTRUCTOR.

//------------------------------------------------------------------------------------------
if(SNG || LNG)
	{//--------just output the potential to file to take a look at it

	if( P->L < 5*a )
		{
		cout << "\n WARNING: the system size you have chosen is less than 5 times the particle size ";
		}

	//-----------------------TEST PLOT TO SEE WHAT THE POTENTIAL LOOKS LIKE --------------
	/*
	clear_charray(cpath, charlength );

	sprintf(cpath, "%sv2%spotentialfull_E0_%lf_a_%d.txt",pathout.c_str(),NGtype.c_str(),P->E0,a);
	v2pot_output.open(cpath);
	for(j=0;j< a;j++)
		{v2pot_output << (j+1) << "\t" << P->v2[j] << endl ;}	
	v2pot_output.close();
		//---------------   TEST PLOT FINISHED HERE    ----------------------------
	*/
	}

	//------------------------------------------------------------------------------------------


(*P).log = flog;

*flog << "\n After ODEdat constructor, Llim is NOW set to " << P->L << endl;

//------------read in parameters from previous iteration--------------

clear_charray(cpath, charlength );
if(should_export_IC)
	{
	sprintf(cpath, "%svoiddat%sBZ%s-initseg_E0-%.2f_muN-%.2f_t_rho_mean_stddev_rhodotnum.txt", pathout.c_str(), NGtype.c_str(), BZ.c_str(),P->E0,P->muN);
	}
	else
	{
	sprintf(cpath, "%svoiddat%sBZ%s-finalseg_E0-%.2f_muN-%.2f_t_rho_mean_stddev_rhodotnum.txt", pathout.c_str(), NGtype.c_str(), BZ.c_str(),P->E0,P->muN);
	}


foutmain = new ofstream(cpath);	

// this is what gets plotted:
// *foutmain << t << "\t" << ((*P).rho)/P->L) << "\t" << P->mean << " \t " << P->std_dev << "\t" << rhodot_num << endl;
//--------------------------------


//-------------------- THE MAIN LINE IN THIS PROGRAM IS HERE: -----------------------------------

			test_result=time_go(P, foutmain);

//------------------------THAT WAS THE MAIN LINE IN THIS PROGRAM---------------------------------


// *flog << "\n upon completion: L=" << L << " rho=" << "; sigma=" << P->std_dev << "; mu="<<P->mean <<"; tf=" << P->t;


cout << "\n upon completion, system density was:" <<  (P->rho)  << " at time" << P->t;
cout << "\n void numerics program complete. exiting successfully.\n";
	
if( ! should_export_IC)
	{
	*maxrhovals  << (P->muN) << " \t " << (P->E0) << " \t " << (P->maxrho/P->L)  << " \t " << (P->rho/P->L) << endl; 
	}
//-------------------------- CLOSE FILES AND CLEAN UP MEMORY ---------------------
(*foutmain).close();
delete foutmain;


delete P;


//---------------- @@@ DOWN TO HERE ---------------------
//--------------------------------------------------------------------------------

(*flog).close();
delete flog;


return 0;
}
//****************=====----- END OF MAIN ---=====*******************


