//---void.h -function definitions for the void numerics calculation.
// ---last updated on  Sun Apr 27 21:18:49 CEST 2014  by  bren  at location  , bren-Desktop

//  changes from  Sun Apr 27 21:18:49 CEST 2014 : implemented criteria to terminate the run when a maximum in the density is encountered -called "should_peakterminate"

//  changes from  Fri Apr 25 17:08:48 CEST 2014 : optimized the process by iteratively cutting off large voids that have gone negative. cutoff is in func_SLNG"

//  changes from  Sun Apr 13 23:34:52 CEST 2014 : modified the script for explicitly-sized particles (i.e. suitable for dimers, trimers, etc., no longer just coarse-grained full-sized Nucl's.)

//  changes from  Tue Nov 26 19:18:48 CET 2013 : ran valgrind to clean up mem. management. Added a destructor to the ODEdata class.
//-------------------------------------------------------------------------------

//-------------this is only the only file that contains explicit differences for the BZ_on_add/removal condition.

 #ifndef __void_h  //---check whether it's been defined already so we don't do it twice. 
 #define __void_h 

#include <fstream>
#include <cstdio>
#include <cstdlib>
   
#include <iostream>  
#include <iomanip>  
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <queue>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <sstream>
#include <bren_lib.h>


using namespace std;


//****************************************************************
const int Np10 = 10; 		// Resolution of time-samples taken. this is probably high enough resolution on


//****************************************************************
ODEdat::ODEdat(const double * rates_times, const  int * sizes, const bool * B, const string pathout_in) //--the constructor.
{
int i; 
//------------------------------------
t0  = rates_times[0];	//---this is always the same: it sets the step-size. 
			// ONLY t is changed based on read-in condition.
t1      = rates_times[1];
rm      = rates_times[2];
muN     = rates_times[3];
E0      = rates_times[4];
t_export = rates_times[5];

L              = sizes[0];	phys_bound = L;	//--- initially.
a              = sizes[1];


shouldplotvoiddist   = B[0];
shouldplotrhos       = B[1];

HNG                  = B[2];
SNG                  = B[3];
LNG                  = B[4];
boltzmann_on_add     = B[5];
boltzmann_on_removal = B[6];
boltzmann_on_uphill  = B[7];
should_check_neg     = B[8];
should_import_IC     = B[9];
should_export_IC     = B[10];
should_peakterminate = B[11];


if( ( int(boltzmann_on_add)+int(boltzmann_on_removal) + int(boltzmann_on_uphill) ) != 1)  
	{
	cout << "\n ERROR: the respective boltzmass conditions is not unique and/or sufficient.\n";
	exit(1);
	}
//------------------------------------------------------

size_v2 =0;

if ( HNG && SNG)
	{
	cout << "\n ERROR: HNG and SNG cannot BOTH be true. \n";
	exit(1);
	}

//------------------------------------

pathout = pathout_in;
plotnum=0;
plotnum_initial=0; 	//--- this will be over-written if import_IC gets called.
attempt=0;

mean 		 = 0.0;
std_dev 	 = 0.0;

//---------------- IMPORT/EXPORT CRITERIA GOES HERE ----------------------

t=0;	//--- this will be overwritten inside import_IC() if we call it.


int    total_obs_vdist;



if(should_import_IC)
	{
	import_IC(); // this will automatically change 'L' and phys_bound and initialize V_IC. Everything else should remain the same.

	total_obs_vdist   = ceil(   Np10* (log10( t1/t) ) )  ;
	}
else
	{//--- no IC to import, then just set the whole system to empty.
	V_IC = new double[L+1];
	for(i=0;i<L;i++)
		{
		V_IC[i]=0.0;
		}
        V_IC[L]=1.0; 
	total_obs_vdist   = ceil(   Np10* (log10( t1/t0) ) )  ;
	}


//--------------------ALLOCATE AND ASSIGN THE 2-BODY INTERACTION -------------------
xcoarse       = new double[L+1];
has_been_neg  = new bool[L+1]; 


has_been_neg[4]=true; // ---just testing my own initialization script.
init_array( has_been_neg , L+1, false );

//--------------coarse graining here:--------------------------------------

rp = gsl_sf_exp(muN); 
r  = rp/rm;		//---keep these two lines constant for now!

if(SNG || LNG )
	{
	//-----this is still valid regardless of when we TAKE the boltzmann factor.
	//-----we're just determining what it is exactly.

	Bzman_v2 = new double[L+1];
	v2       = new double[L+1];

	init_array(v2      , L+1, 0.0);
	init_array(Bzman_v2, L+1, 0.0);
	init_array(xcoarse , L+1, 0.0);	//---initialization functions are defined in "bren_lib.h"

	init_array( v2 , a, 0.0);

	if(SNG)
		{
		VNN_SNG_calc_smallp(v2, L, a, E0); 
		}
	else if(LNG)
		{
		VNN_LNG_calc_smallp(v2, a, E0); 
		}

	//---- got the potential: from here, SNG and LNG are treated the same 
	//------(they just have a different v2 potential from here).

	for(i=0;i<=L;i++)
		{
		Bzman_v2[i] = gsl_sf_exp(-v2[i]);
		}
	}
else if (HNG)
	{
	Bzman_v2 = NULL;	//----- code it this way to act as a warning in case v2 
	v2       = NULL;	//----- is ever accessed during an HNG run (it shouldn't be).

	/*-------here is where we hard code the HNG in the same way as SNG.------
	for(i=0;i<=L;i++)
		{
		if(i<a)
			{
			Bzman_v2[i] = 1.0;
			}		
		else if(i>=a)
			{
			Bzman_v2[i] = 0.0;
			}		
		}
	-------------------------------------------------------------------------*/
	}

//-------------------- logarithmic spacing of observation points -------------------------

printq        = new bool[total_obs_vdist];
tpoints_vdist = new double[total_obs_vdist];


space_tpoints_logarithmic( t, t1, Np10 , total_obs_vdist, tpoints_vdist);

init_array( printq,  total_obs_vdist , false );

/*** --- CANT BE BOTHERED BEING FANCY ABOUT THIS, JUST MAKE THEM EVENLY SPACED ------

bool found=false;
for(i=0;i<Np10;i++)
	{
	if (t <= dummy_t_points[i])
		{
		found = true;
		break;
		}
	}
if ( found == false)
	{
	cout  << "\n ERROR: cannot find t0>dummyt0[i]. Exiting. \n ";
	exit(1);
	}

int dummy_num_cutoff = i;

total_obs_vdist = dummy_num_t_points-i;	//--- cut off the leading points that come before t0


for(i=0;i<total_obs_vdist;i++)
	{
 	printq[i]=false;
	tpoints_vdist[i] = dummy_t_points [ i + dummy_num_cutoff] ;
	if ( i + dummy_num_cutoff >= dummy_num_t_points )
		{ 
		cout << "\n ERROR: reading off array end dummy_t_points."; exit(1); 
		}
	}
 -------------------------------------------------------------------------*/


//----  @@@ no longer simply space from t0 like this: 
//----  space_tpoints_logarithmic(t0, t1, Np10 , total_obs_vdist, tpoints_vdist );
//----------------------------------------------------------------------------------------


rho     = 0.0;
maxrho  = 0.0;
rhostar = 0.0;
C1      = 0.0;
C2      = 0.0;
//! rhodot  = 0.0;
}

//************************************************************************END OF ODEdat constructor
ODEdat::~ODEdat() //destructor 
{
delete [] xcoarse;

if(SNG || LNG )
	{
	delete [] Bzman_v2;
	delete [] v2;
	}

delete [] printq;
delete [] tpoints_vdist;
delete [] has_been_neg;
delete [] V_IC;		//---this is ALWAYS created and must always be destroyed. It just happens to be all zeros if there's no import.

}

//*************************************************************************

//{}
//int jac (double t, const double c[], double *dfdy, double dfdt[], void *params)
//*************************************************************************
string convertDouble(double number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

//****************************************************************
ODEdat& ODEdat::operator= (const ODEdat& param)
{
v2          	= param.v2;		//----- v2[x]       = the energetic potential of two particles with 'x' empty sites between them.
Bzman_v2	= param.Bzman_v2;	//----- Bzman_v2[x] = exp(-v2[x]), the boltzmann factor thereof.
size_v2		= param.size_v2;		// the number of points in the v2 potential (beyond which it is assumed v2=0)

rho		= param.rho;		//---the average density.
rhostar		= param.rhostar;	//---the normalization factor.
maxrho		= param.maxrho;		//---the derivative
empty_space	= param.empty_space;
C1		= param.C1;
C2		= param.C2;

rm		=param.rm; // = r-'minus' -the off rate.
rp		=param.rp; // = r-'plus'  -the on rate.
r		=param.r;  // = r+/r-

a		=param.a; //---the finite size of the particles.
L		=param.L; //---the size of the system /max void.
t		=param.t;
t1		=param.t1;

printq		=param.printq;		//----the queue of plots to make.
total_obs_vdist	=param.total_obs_vdist;		//----the number of plots of the void distribution we should make.
shouldplotvoiddist	=param.shouldplotvoiddist; //---whether or not we should plot the void distribution at various times.
shouldplotrhos		=param.shouldplotrhos;
plotnum			=param.plotnum; //----the current index of the void-distribution plot.

cpath			=param.cpath;   // [charlength];

attempt			=param.attempt;
pathout			=param.pathout;

HNG			=param.HNG;
SNG			=param.SNG;

log			=param.log;
step			=param.step;

return *this;

}

//****************************************************************

int ODEdat::get_mean_stddev( const double * V )
{
int i=0,j=0;
bool result;
double temp=0.0;

mean 		 = 0.0;
std_dev 	 = 0.0;
double Voidsum   = 0.0;
double checknorm = 0.0;

for(i=0;i<=L;i++)
	{
	Voidsum += V[i]; //---normalization factor for Void distances;
	}

for(i=0;i<=L;i++)
	{
	mean += i*(V[i]/Voidsum);  //--- i is the void size; V[i]/Voidsum is the probability.
	checknorm += V[i]/Voidsum; //--- normalization factor for Void distances;
	}

for(i=0;i<=L;i++)
	{
	std_dev += ( (V[i]/Voidsum)* (double(i)-mean)*(double(i)-mean) ); //---calculate the std.dev.
	}

std_dev=sqrt(std_dev);

if(fabs(checknorm-1) > 0.0000001)
	{
	cout << "\n ERROR in get_mean_std_dev. checknorm is not checking out.\n\n";
	exit(1);
	}

return result;
}
//****************************************************************

bool ODEdat::printtime(void)
{
int i=0,j=0;
bool result;
double temp=0.0;

while(printq[i] == true)	// printq starts out all false, the first one will be set true only
	{
	i++;		// when t > (1/total_obs)*t1, etc.
	if(i >=total_obs_vdist)
		{
		return false;
		}
	}

temp = tpoints_vdist[i];

if(t > temp)
	{
	printq[i]=true;
	result	=true;
	}
else	
	result = false;


//!----fix this later if you want wf's---
// result = false; //---ALWAYS!
return result;
}

//*********************************************************************************
int ODEdat::printout_voiddist(const double * V)
{
int NV;
int res;

int i, j;
int  charlength=400; //--the number of characters in the string for our path.
char cpath[charlength];
clear_charray(cpath, charlength );


//----------------------------------------------------------------
if( plotnum < 10)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_0000%d.txt", pathout.c_str(),plotnum+plotnum_initial);
	}
else if(plotnum < 100)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_000%d.txt", pathout.c_str(),plotnum+plotnum_initial);
	}
else if(plotnum < 1000)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_00%d.txt", pathout.c_str(),plotnum+plotnum_initial);
	}
else if(plotnum < 10000)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_0%d.txt", pathout.c_str(),plotnum+plotnum_initial);
	}
else 
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_%d.txt", pathout.c_str(),plotnum+plotnum_initial);
	}

ofstream fvout(cpath);

/*-----------NO NEED TO BOTHER WITH THIS SHIT ANYMORE----------------
//-------get Redner's formula:
double A     = coverage*coverage*(CGF*L) /( kHNG_original * (coverage + kHNG_original*(1.0-coverage )  ) );
double alpha = gsl_sf_log( 1 + coverage/(kHNG_original*(1.0-coverage) )  );
double A_CG     = coverage*coverage*(L) /( a * (coverage + a*(1.0-coverage )  ) );
double alpha_CG = gsl_sf_log( 1 + coverage/(a*(1.0-coverage) )  );
*/

//--------------------------PLOT THE VOID DISTRIBUTION AS A FUNCTION OF X ------------------------
for(i=0;i<L;i++)
	{
	fvout << (i+1) << " \t " << V[i]/(L) << " \t " << V[i]/rho << " \t " << has_been_neg[i] << endl;
	}
//---------------NOW THE TIME ASSOCIATED WITH THIS PARTICULAR plotnum iteration------------------

sprintf(cpath, "%splotnum_t_rho_coverage.txt", pathout.c_str());
ofstream pn;

if(plotnum==0)
	pn.open(cpath);
else
	pn.open(cpath,fstream::app);

pn << plotnum << " \t " << t << " \t " << (rho/L) << " \t " << coverage << endl;

//------------------------------------------------------------------------------------------------

pn.close();
fvout.close();

return 1;
}
//*********************************************************************************

/* decided not to use this function anymore, just always start from empty. 

int ODEdat::initialize_V(const int init, double * V)
{
// V[NV-1]=1; //---start from empty, meaning only the largest possible void has density one.
// V[init]=(float(P->L)/float(init+P->k)); //---start from full
// V[init]=1;
// V[P->L -2*P->k -1]=1;
// V[init]=1;

int Qi = floor(L/(init+k));//--the number of particle, m-void combinations that can fit.

int remainder = L%(init + k);
int space_leftover=0;


if (remainder == 0)	
	{
	V[init] = floor(L/(init + k)); 	//----easy peasy: system size works out to an 
					//----integer multiple of the void-particle combination

	}	
else
	{				//---slightly more complicated
	
//  0                                                             L
//  |-------------------------------------------------------------|
//  |+++|<---->|+++|<---->|+++|<---->|+++|<---->|+++|<---->...
//            1|         2|         3|         4|         5|    Qi=5
//  |                                               |<----------->|  				       
//  Qi is the number of particles. here, Qi=5,        ^but this is now "space leftover"

if(Qi >=1 ) //---then just subtract off the last one
	{
	V[init]           = Qi-1; //the number of voids of our interest is the particle number minus one.
	space_leftover    = L-((Qi-1)*(init+k) +k);
	V[space_leftover] = 1.0;
	}
else
	{
	V[L]=1;//----completely empty system.
	}

 	*log << "\n Site distribution works out unevenly. residual placed in voids of size " << space_leftover << endl;

	}

 *log << "for init, we have set V[" << init << "] = " << V[init] << ", V[" << space_leftover << "] = " << V[space_leftover] << endl;

 *log << "\n ------------------------------------------";
return 1;
}
!!!!!*/


//*********************************************************************************
int func_HNG (double t, const double V[], double f[], void *params)
{
//--- THIS FUNCTION IS WHERE THE BULK OF THE COMPUTATION TAKES PLACE.----------
//--- THE ARRAY OF DERIVATIVES f[] IS CALCULATED USING THE MATRIX ELEMENTS V[].

ODEdat*  P = (ODEdat *)params;
//P->V = V;

int j,m; //-----counter ints
double C=0.0;
int  kn = P->a; // ----the particle size.
int  L  = P->L; // ----the system size.
//---------------------------------------------------------
P->attempt +=1;		// keep track of the number of times we've "attempted to time-step here"
P->t = t;		// the current time.

double rm = P->rm;
double rp = P->rp;	// rp/m = 'plus/minus' --the on off rates.

int  phys_bound = P->phys_bound;

//--------------------  initialize  ----------------------

init_array( f,  L+1, 0.0 );



// (*P).get_rho_anal(V);  //---returns exactly the number of particles in the system -WE DON'T WORRY ABOUT COARSE-GRAINING HERE!
(*P).get_rhostar(V); 
(*P).get_empty_space(V);
 
//----make sure to update rho to fit the numerics in the convolution term (number 4)

//============= here handle finite-size special cases =================================================

if( ! P->has_been_neg[L] )
	{
	//----completely empty system  		0|---------------------------------|L	
	f[L] = -1*rp*L*V[L] + rm*V[L-kn] ;	//--only terms 3,4 survive here.
	}

if( ! P->has_been_neg[L-kn] )
	{
	//----single particle in system         0|------------|<-kn->|--------------|L	
	f[L-kn]  = -1*rm*V[L-kn] + rp*L*V[L];     //---term 1 only has a single rm-removal rate., 
					  	  // term 2 has full range of 'L' binding sites.

	if ( (L-2*kn+1) >=1)
		{//---term 3 : same as canonical case.
		f[L-kn] -= rp*(L-2*kn+1)*V[L-kn];    
		}

	//---- term (4) ---creation of void 'm' from desorption 
	//---- of neighbour with next void adding up to m.   

	for(j=0; j<= (L-2*kn);j++)	 
		{ 
		if( !isnan(1.0/(P->rhostar)) && !isinf(1.0/(P->rhostar)) ) //---- float-error catch (in case rho=0)
			{
			f[L-kn] += (rm/(P->rhostar))*V[j]*V[L-2*kn-j];
			}
		}

	}

//==== done handling special finite-size cases; the rest are canonical ===================
//--- optimized out: for(m=0; m<= (L-2*kn) ;m++)

for(m=0; ( (m<= (L-2*kn)) &&(m<=phys_bound) ) ;m++)
	{
	//----term(1) -particles on either side leaving.
	f[m]=-2*rm*V[m]; 
	//----term (2) ---creation of void 'm' from a larger one.  
	for(j=m+kn; j<=(L-kn) ;j++)	
		{ 
		f[m] += 2*rp*V[j]; //---term 2
		}

	if(m>=kn)
		{
		//---term(3) disappearance of V[m] due to adsorption in its interior.	
		f[m] += -rp*(m-kn+1)*V[m];

		//----term (4) ---creation of void 'm' from desorption of neighbour with next void adding up to m.   
		for(j=0; j<= (m-kn);j++)	 
			{ 
			if( !isnan(1/(P->rhostar)) && !isinf(1/(P->rhostar)) ) //---- in case rho=0
				{f[m] += (rm/(P->rhostar))*V[j]*V[m-kn-j];}
			}
		}

       }//--- end the if/else special conditions of m.

return GSL_SUCCESS;
}
//*********************************************************************************
int func_SLNG (double t, const double V[], double f[], void *params)
{
//--- THIS FUNCTION IS WHERE THE BULK OF THE COMPUTATION TAKES PLACE.----------
//--- THE ARRAY OF DERIVATIVES f[] IS CALCULATED USING THE MATRIX ELEMENTS V[].

//--- HERE WE CONTAIN EXPLICIT DIFFERENCES FOR THE BZ_on_add/removal CONDITION.


ODEdat*  P = (ODEdat *)params;
//P->V = V;

int i=0,j=0,m=0; //-----counter ints
double C=0.0;
int  L = P->L; // ----the system size.
int  phys_bound = P->phys_bound;

//---------------------------------------------------------
P->attempt +=1;		// keep track of the number of times we've "attempted to time-step here"
P->t = t;		// the current time.

double rm = P->rm;
double rp = P->rp;	// rp/m = 'plus/minus' --the on off rates.

//--------------------  initialize  ----------------------

init_array( f,  L+1, 0.0 );

// (*P).get_rho_anal(V);  //---returns exactly the number of particles in the system
(*P).get_rhostar(V); 
(*P).get_empty_space(V);
 
//----make sure to update rho to fit the numerics in the convolution term (number 4)

if( P->boltzmann_on_add)
	{
	//======= here handle finite-size special cases =========================

	if( ! P->has_been_neg[L] )
		{
		//----completely empty system  		0|---------------------------------|L	
		f[L]= -1*rp*L*V[L]*P->Bzman_v2[L-1] + rm*V[L-1] ; //--only terms 3,4 survive here.
		}

	if( ! P->has_been_neg[L-1] )
		{
		//----single particle in system         0|------------|<-k->|--------------|L	
		f[L-1]  = -1*rm*V[L-1] + 1*rp*L*V[L]*P->Bzman_v2[L-1];     
		//---term 1 only has a single rm-removal rate., 
					// term 2 has full range of 'L' binding sites.
		for(i=0;i<(L-1);i++)
			{
			f[L-1] -= rp * V[L-1] * ( (P->Bzman_v2[i]*P->Bzman_v2[L-1-1-i] )/ (P->Bzman_v2[L-1]) )   ;    
			//---term 3 : same as canonical case.
			}

		//----term (4) ---creation of void 'm' from desorption
		//----  of neighbour with next void adding up to m.   

		for(j=0; j<= (L-1-1);j++)	 
			{ 
			if( !isnan(1.0/(P->rhostar)) && !isinf(1.0/(P->rhostar)) ) 
				{ //---- THIS IS A float-error catch (in case rho=0)
				f[L-1] += (rm/(P->rhostar))*V[j]*V[L-1-j];
				}
			}
		}		
	//---- done handling special finite-size cases; the rest are canonical ==========
	for(m=0; ((m < (L-1)) &&(m<=phys_bound) ) ;m++)
		{
		//----term(1) -particles on either side leaving.
		f[m]=-2*rm*V[m]; 

		//----term (2) ---creation of void 'm' from a larger one.  
		for(i=m+1; (i<=(L-1) && i <= phys_bound ); i++)	// --- V'>phys_bound are all 0.
			{ //---term 2 
			f[m] += 2*rp*V[i] * ( (P->Bzman_v2[m]*P->Bzman_v2[i-1-m] )/ (P->Bzman_v2[i]) ); 
			}

		//---term(3) disappearance of V[m] due to adsorption in its interior.	
		for(i=0; i<=(m-1) ;i++)	
			{ //---term 3----
			f[m] -= rp*V[m] * ( (P->Bzman_v2[i]*P->Bzman_v2[m-1-i] )/ (P->Bzman_v2[m]) ); 
			}

		//----term (4) ---creation of void 'm' from desorption of neighbour with next void adding up to m.   
		for(j=0; j<= (m-1);j++)	 
			{ 
			if( !isnan(1/(P->rhostar)) && !isinf(1/(P->rhostar)) ) //---- in case rho=0
				{
				f[m] += (rm/(P->rhostar))*V[j]*V[m-1-j];
				}
			}
       		}//--- end the if/else special conditions of m (i.e. done with the canonical cases.)

	}//-----close the "if(boltzmann_on_add)" condition

else if( P->boltzmann_on_removal)//-------------------bolztmann factors on removal--------------
	{

	//	@@@ STILL NEED TO DO THE SAME THING HERE

	// N.B. the Bzman array is always exp(-v2(x)), and therefore always <=1.
	// -in order to speed up removal, we have to DIVIDE rates by it.

	if( ! P->has_been_neg[L] )
		{
		//---- L ------------------------------------------	
		f[L]= -1*rp*L*V[L] + rm*V[L-1]/(P->Bzman_v2[L-1]) ; //--only terms 3,4 survive here.
		}

	if( ! P->has_been_neg[L-1] )
		{
		//---- L-1 ----------------------------------------
		f[L-1]  = -1*rm*V[L-1]/(P->Bzman_v2[L-1]) + 1*rp*L*V[L];   //--- this includes both (1) loss due to 
								   //--- removal and (2) creation of [L-1] 
								   //--- from larger->addition.

		f[L-1] -= rp*(L-1) * V[L-1];    //--- term 3 : same as canonical case 
					 //--- equal binding possibilities everywhere.

		for(j=0; j<= (L-1-1);j++)	//--- term 4 convolution 
			{ 			//
			if( !isnan(1.0/(P->rhostar)) && !isinf(1.0/(P->rhostar)) ) 
				{ //---- THIS IS A float-error catch (in case rho=0)
				f[L-1] += (rm/(P->rhostar))*V[j]*V[L-1-1-j] * ( (P->Bzman_v2[L-1])/ (P->Bzman_v2[j]*P->Bzman_v2[L-1-1-j] ) )    ;
				}
			}
		}
	
	
	//---- done handling special finite-size cases; the rest are canonical ==========

	//------- optimized out: for(m=0; m< (L-1) ;m++)
	for(m=0; ((m < (L-1)) &&(m<=phys_bound) ) ;m++)
		{
		//----term(1) -particles on either side leaving.
		for(i=0; i<=(L-2-m) ;i++)	
			{
			f[m] -= 2*rm*V[m] * (V[i]/P->rhostar) * (1/(P->Bzman_v2[m]*P->Bzman_v2[i] ))*P->Bzman_v2[i+m+1] ; 
			}
		//----term (2) ---creation of void 'm' from a larger one.  
		// @@@ --- optimized out: for(i=m+1; i<=(L-1) ;i++)	

		for(i=m+1; (i<=(L-1) && i <= phys_bound ); i++)	// --- V'>phys_bound are all 0.

			{ //---term 2 
			f[m]+=2*rp*V[i]; 
			}
		
		//---term(3) disappearance of V[m] due to adsorption in its interior.	

		f[m] -= rp*V[m] * m; 

		//----term (4) ---creation of void 'm' from desorption of neighbour with next void adding up to m.   
		for(i=0; i<= (m-1);i++)	 
			{ 
			if( !isnan(1/(P->rhostar)) && !isinf(1/(P->rhostar)) ) //---- in case rho=0
				{
				f[m] += (rm/(P->rhostar))*V[i]*V[m-1-i]* (P->Bzman_v2[m] /(P->Bzman_v2[i]*P->Bzman_v2[m-1-i])) ;
				}
			}

       		}//--- end the if/else special conditions of m (i.e. done with the canonical cases.)
	}


return GSL_SUCCESS;
}

//*********************************************************************************
double ODEdat::get_rhostar(const double * V)
{//------from page 63 in notebook (the normalization factor that should impose conservation criteria)
int j,n,m;
double result= 0.0;
double denom=0.0;
double num=0.0;

double inner_sum=0.0;

if(HNG)
	{
	//------------ GET THE HNG NUMERATOR -------
	for (n=0;n<=(L-3*a);n++)
		{
		inner_sum = 0.0;

		for (j=0;j<=n;j++)
			{
			inner_sum += (V[j]*V[n-j]);
			}

		num += (double(n)+2.0*double(a))*inner_sum;
		}
	n=L-2*a;//=1
	inner_sum = 0.0;

	for (j=0;j<=n;j++)
		{
		inner_sum += (V[j]*V[n-j]);
		}

	num += (double(n)+2.0*double(a))*inner_sum;

	//------------GET THE HNG DENOMINATOR ----------
	for (m=0;m<=(L-2*a);m++)
		{
		denom +=2*(a+m)*V[m];
		}
	//------------------------------------------
	}
else if(SNG || LNG)
	{
	//------------ GET THE SNG NUMERATOR -------
	for (n=0;n<=(L-1);n++)
		{
		inner_sum = 0.0;

		for (j=0;j<n;j++)
			{
			inner_sum += (V[j]*V[n-1-j]);
			}

		num += (double(n)+1.0)*inner_sum;
		}

	//------------GET THE SNG DENOMINATOR ----------
	for (m=0;m<=(L-2);m++)
		{
		denom +=2*(1+m)*V[m];
		}
	//------------------------------------------
	}



if( !isnan(num/denom) && !isinf(num/denom) ) //---- float-error catch (in case rhostar=0)
	result  = num/denom;
else 
	result = 1;

rhostar = result;
return result;
}
//*********************************************************************************

/*
double ODEdat::get_rhodot_anal(const double * V)
{//------from page 28 in notebook (analogous to eq. 7.60 in Redner.)
int j;
double result= -1.0*rho;

if(HNG)
	{
	for(j=k; j<L; j++)
		result += rp*(j-k+1)*(V[j]);

	rhodot=result;
	}
else if (SNG)
	{
	cout << "\n ERROR: haven't resolved this case yet. exiting.\n";
	exit(1);
	}


return result;
}
*/
//*********************************************************************************
double ODEdat::get_rho_anal(const double * V)
{
int i;
double result=0.0;
for(i=0;i<L;i++)
	result += V[i];

rho=result;
if(HNG)
	{
	coverage = rho*a/L;
	}
else if(SNG || LNG )
	{
	coverage = rho*a/L; //--implement the "coverage from DNA perspective later.
	}
else
	{
	cout << "\n ERROR: undefined NG type\n";
	*log << "\n ERROR: undefined NG type\n";
	(*log).close();
	exit(1);
	}

return result;
}
//*********************************************************************************
double ODEdat::get_C(const double * V)//---the conserved quantity. This should always add up to zero.
{
int m,n,j;
double result=0.0;
C1=0.0; C2=0.0;
double inner_sum=0.0;


//----------------------------------------------
for(m=0;m<=(L-2*a);m++)
	{
	C1 += (-2*rm)*(a+m)*V[m];
	}
//----------------------------------------------
for(n=0;n<=(L-3*a);n++)
	{
	inner_sum=0.0;
	
	for(j=0;j<=n;j++)
		{
		inner_sum +=V[j]*V[n-j];		
		}
	C2+= (n+2*a)*rm*inner_sum;
	}
//---------now add on the last one for L-2a-----------
n=L-2*a;
inner_sum = 0.0;

for (j=0;j<=n;j++)
	{
	inner_sum += (V[j]*V[n-j]);
	}
C2+= (n+2*a)*rm*inner_sum;
//----------------------------------------------

result = C1+C2/rhostar;

return result;
}

//*********************************************************************************
double ODEdat::get_empty_space(const double * V)
{ //----this is the same for both cases, HNG or SNG. We define the footprint length outside of the dyad to be "empty"
int i;
double result=0.0;
for(i=0;i<=L;i++)
	result += i*V[i];

empty_space=result;

result=result;
}

#endif //--this ends the clause as to whether or not this symbol (i.e. this file) has been defined already.

//*********************************************************************************
int  ODEdat::import_IC(void)
{
int i=0,x=0;
int  charlength=400; //--the number of characters in the string for our path.
char cpath[charlength];

int    a_imported=0, L_imported=0;
double E0_imported  =0.0;
double muN_imported =0.0;
double t_imported   =0.0;
double rho_imported =0.0;
double dummy=0.0;

//----------------------------------------------------------------
clear_charray(cpath, charlength );
sprintf(cpath, "./IC_import/import_t0_rho0_muN-%.2f_eps-%.2f_a-%d.txt", muN, E0, a);

ifstream t0_rho0_in(cpath);
if( t0_rho0_in.fail() )
	{
	cout << "ERROR: could not open parameter IC file named : ";
	cout << cpath << "\n Now exiting.\n";
	exit(1);
	}

t0_rho0_in >> t_imported >> rho_imported >> dummy >> L_imported  ;
t0_rho0_in.close();

double  V_temp[L_imported+1] ;
int     has_been_neg_temp[L_imported+1]; 

init_array( V_temp,            L_imported+1, 0.0 );
init_array( has_been_neg_temp, L_imported+1, false );


//----------------------------------------------------------------
clear_charray(cpath, charlength );
sprintf(cpath, "./IC_import/import_V0dist_muN-%.2f_E0-%.2f_a-%d.txt", muN, E0 , a );

ifstream vdist_in(cpath);
if (vdist_in.fail() )
	{
	cout << "\n ERROR: failed to access IC import file.\n" ;
	exit(1);
	}

x=0;
double Norm=0.0;
double Norm_inc=0.0;

int has_been_neg_check;
double Vdenscheck=0.0, Norm_inc_check=0.0;

for(i=0;i<=L_imported;i++)
	{
	vdist_in >> x >> V_temp[i] >> Norm_inc >> has_been_neg_temp[i];
	Norm += Norm_inc;
	if( has_been_neg_temp[i] || x!=i+1)
		{
		goto error_reading_input_Vdist;
		}
	}

vdist_in >> x >> Vdenscheck >> Norm_inc_check >> has_been_neg_check;
vdist_in.close();

if( !has_been_neg_check ||  Norm_inc_check != 0.0 || Vdenscheck!= 0.0 || Norm-1.0 >0.00001 || Norm-1.0 < - 0.00001 )
	{
	error_reading_input_Vdist:
	cout << "\n ERROR: SOMETHING went wrong when reading in the Vdistribution file in import_IC. Exiting. \n";
	exit(1);
	}

//--- L_imported is now the last index which is ruled physical.
//--- THIS IS THE SIZE OF THE NEW SYSTEM, AND WE MUST ALLOCATE ONE MORE FOR THE ZERO PLACE.

//----------------------------------------

L           = L_imported;
phys_bound  = L_imported;

V_IC = new double[L+1];
for(i=0;i<=L;i++)
	{
	V_IC[i] = L * V_temp[i];
	}
//----------------------------------------------------------------
t   = t_imported;
rho = rho_imported*double(L);

}

//*********************************************************************************
int  ODEdat::export_IC(double * V)
{
int i;
int  charlength=400; //--the number of characters in the string for our path.
char cpath[charlength];
double dummy = 0.0;
double 	rho_t	    = get_rho_anal(V);
//----------------------------------------------------------------
clear_charray(cpath, charlength );
sprintf(cpath, "%simport_t0_rho0_muN-%.2f_eps-%.2f_a-%d.txt", pathout.c_str(), muN, E0, a);
ofstream params_out(cpath);


params_out << std::setprecision(10) << t << " \t " << rho/double(L) << " \t " << dummy << " \t " << phys_bound  << endl;
params_out.close();
//----------------------------------------------------------------
clear_charray(cpath, charlength );
sprintf(cpath, "%simport_V0dist_muN-%.2f_E0-%.2f_a-%d.txt", pathout.c_str(), muN, E0, a);

ofstream vdist_out(cpath);

for(i=0;i<=L;i++)
	{
	vdist_out << (i+1) << " \t " << V[i]/L << " \t " << V[i]/(rho_t) << " \t " << has_been_neg[i] << endl;
	}
vdist_out.close();
//----------------------------------------------------------------

}
