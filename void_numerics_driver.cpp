/*--- VOID NUMERICS DRIVER SCRIPT -TESTS THE EQUATIONS DERIVED FOR SNG IN ANALOGY WITH REDNER -----
// ---last updated on  Tue Nov 26 19:18:48 CET 2013  by  Brendan.Osberg  at location  th-ws-e537

//  changes from  Tue Nov 26 19:18:48 CET 2013 : ran valgrind to clean up mem. management. Added a destructor to the ODEdata class.
//-----------------------------------------------------------------------------------------------*/


#include <fstream>
#include <cmath>
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
int     L=0, w=0, j=0, test_result=0; 
//---the size of the system, # of plots to make, size of particle, size of footprint, counting index, dummy.
double  t0, t1, rm, E0, CGF, kHNG_CG, kHNG_exact_test, irho_target , muN_input;
bool    shouldplotvoiddist, shouldplotrhos, HNG, SNG, LNG;
bool    boltzmann_on_uphill, boltzmann_on_add, boltzmann_on_removal;

string  pathout;

double muN_original=0.0; // chemical potential used to determine rp;

//--'command line parameters'------------------------------------------------------------------------
if (argc != 3)
	{
	cout << "\n error: expecting 2 command-line argument (1) taskID and (2) NG type. Instead we got " << argc << " arguments which were:\n";
	for (j=1;j<argc;j++)
		{
		cout << argv[j] << ",  ";
		}
	cout << "\n exiting\n\n";
	exit(1);
	}
int TASKID    = atoi(argv[1]); //---atof is for floats and doubles.
string NGtype = argv[2];

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

fin  >> muN_input >> irho_target; 	// chemical potential without coarse-graining. *MAYBE* the irho_target parameter is used.
fin  >> kHNG_exact_test >> kHNG_CG;	// size of the particle before and after coarse-graining. should be usually 147
fin  >> w >> E0;	// size of the footprint

fin  >> shouldplotvoiddist >> shouldplotrhos;		// boolean should we print the void dist profiles.
fin  >> BZflag;
fin  >> parity_check;
fin  >> pathout;
fin.close();

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
	rm=0.0;		// rate off ALWAYS =1 
	}
else
	{
	rm=1.0;		// rate off ALWAYS =1 
	}




//---------------------------------------------------------


CGF=kHNG_exact/kHNG_CG; //------DO THE COARSE-GRAINING CALCULATION ON OUR OWN HERE.
string CGF_getmu_feedin;

// /*--------@@@  

irho_target = irho_target + (double(TASKID-1))*5.0;

if(kHNG_CG == kHNG_exact)
	{
	CGF_getmu_feedin ="1";	
	}
else
	{
	clear_charray(cpath, charlength );

	sprintf(cpath, "%do%d",kHNG_exact,int(kHNG_CG));
	CGF_getmu_feedin = cpath;
	}

cout << "\n the CGF string is :" << CGF_getmu_feedin << endl;
muN_original = interpolate_mu_from_rhoi(irho_target,w,E0,kHNG_exact,NGtype, CGF_getmu_feedin );

// -----@@@ **/

/* ! @@@
muN_original = muN_input + 1.0*(TASKID-1); 
@@@ ! */

cout << "\n choice of target irho=" << irho_target << ", leads to selection of muN=" << muN_original << endl;

double muN_CG = muN_original + gsl_sf_log(CGF); //---scale the effective chemical potential according to the coarse-graining.

//---this is done for both SNG and HNG cases.


int init=L;
if (kHNG_exact_test != kHNG_exact )
	{
	cout << "\n ERROR in comparison of expected kHNG_exact. \n"; 
	exit(1);
	}

if (parity_check != 88855888)
	{
	cout << "\n ERROR in ordering of input parameters. exiting. \n"; 
	exit(1);
	}
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
	if ( fabs (kHNG_CG -round(kHNG_CG) ) > 0.0000001) //----make sure it's an integer
		{
		cout << "\n ERROR: k value is not an integer. Exiting.\n\n";
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

	if( ( (SNG || LNG) && L < 5*(2*w+1)/CGF)  || (HNG && L < 5*kHNG_CG ) )
		{
		cout << "\n WARNING: system size is less than 5 particles. exiting.\n";
//!		cout << " the interaction range of the coarse-grained particles. too small. Exiting.\n";
//!		exit(1);
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

if(HNG)
	{
	*flog << "\n new k (effective) is: " << kHNG_CG;
	*flog << "\n L/k ratio  =" << (L/kHNG_CG) << endl;
	}
else if(SNG || LNG )
	{
	*flog << "\n new footprint p is "   << (2*float(w)/CGF);
	*flog << "\n L/footprint ratio is " << L/(2*float(w)/CGF);
	}
//-----------------------------------------------------------------------
*flog << "\n coarse-graining is set to : " << CGF;
*flog << "\n based on mu=" << muN_original;
*flog << "\n chosen for irho_target =" << irho_target;
*flog << "\n system size, L=" << L;
*flog << "\n NGtype = " << NGtype;
*flog << "\n t0 = " << t0 ;
*flog << "\n t1 = " << t1 ;
*flog << "\n kHNG_exact = " << kHNG_exact;
*flog << "\n kHNG_CG    = " << kHNG_CG;
*flog << "\n w = " << w << ", E0 = " << E0;
*flog << "\n shouldplotvoiddist = " << shouldplotvoiddist << ", shouldplotrhos = " << shouldplotrhos;
*flog << "\n BZflag =" << BZflag << endl;
*flog << "\n choice of target irho=" << irho_target << ", leads to selection of muN=" << muN_original << endl;
*flog << "\n CGF = " << CGF << endl;

//-----------------------------------------------------------------------

double rates_times[6];
rates_times[0] = t0;
rates_times[1] = t1;
rates_times[2] = rm;
rates_times[3] = muN_original;
rates_times[4] = E0;
rates_times[5] = CGF;	//----coarse-graining factor.


int sizes[3];
sizes[0] = L;
sizes[1] = kHNG_CG;	//---the finite size of the particles in the coarse-grained system.
sizes[2] = w;


bool B[8];
B[0] = shouldplotvoiddist;
B[1] = shouldplotrhos;
B[2] = HNG;
B[3] = SNG;
B[4] = LNG;
B[5] = boltzmann_on_add;
B[6] = boltzmann_on_removal;
B[7] = boltzmann_on_uphill;

//------------------------------------------------------------------------

ODEdat* P;   //----'P' is simply our catch-all data structure.
ofstream *foutmain;
ofstream test;


//----this is what gets plotted:
//*foutmain << t << "\t" << (((*P).rho)/(P->L*P->CGF)) << "\t" << P->mean << " \t " << P->std_dev << "\t" << rhodot_num << endl;
//--------------------------------

P=new ODEdat(rates_times, sizes, B, pathout);   // ---THE 2-BODY INTERACTION IS CALCULATED 
						// ---AND ASSIGNED IN THE CONSTRUCTOR.


if(SNG || LNG)
	{//--------just output the potential to file to take a look at it

	if( P->L < 5*(2*w+1)/CGF )
		{
		cout << "\n WARNING: the system size you have chosen is less than 5 times ";
		}


	//-----------------------TEST PLOT TO SEE WHAT THE POTENTIAL LOOKS LIKE --------------
	clear_charray(cpath, charlength );

	sprintf(cpath, "%sv2%spotentialfull_E0_%lf_w_%d.txt",pathout.c_str(),NGtype.c_str(),E0,w);
	test.open(cpath);
	for(j=0;j<(2*w+1);j++)
		{test << (j+1) << "\t" << P->v2_uncoarsened[j] << endl ;}	
	test.close();
	clear_charray(cpath, charlength );

	sprintf(cpath, "%sv2%spotential_coarsening_by_%lf_E0_%lf_w_%d.txt",pathout.c_str(),NGtype.c_str(),CGF,E0,w);
	test.open(cpath);

	for(j=0;j<=L;j++)
		{test << P->xcoarse[j]  << "\t" << P->v2[j] << endl;}	
	test.close();
	//------------------------------TEST PLOT FINISHED HERE --------------------------------*/
	}

clear_charray(cpath, charlength );
sprintf(cpath, "%svoiddat%sBZ%s_t_rho_mean_stddev_rhodotnum.txt", pathout.c_str(), NGtype.c_str(), BZ.c_str());

foutmain = new ofstream(cpath);	


(*P).log = flog;

//-



//------------read in parameters from previous iteration--------------

test_result=time_go(P, foutmain);

*flog << "\n upon completion: L=" << L << " rho=" <<  (((*P).rho)/(P->L*P->CGF))  << "; sigma=" << P->std_dev << "; mu="<<P->mean <<"; tf=" << P->t;
(*flog).close();


cout << "\n upon completion, system density was:" <<  (((*P).rho)/(P->L*P->CGF))  << " at time" << P->t;
cout << "\n void numerics program complete. exiting successfully.\n";

//-------------------------- CLOSE FILES AND CLEAN UP MEMORY ---------------------
(*foutmain).close();
delete foutmain;

(*flog).close();
delete flog;

delete P;
//--------------------------------------------------------------------------------


return 0;
}
//****************=====----- END OF MAIN ---=====*******************


