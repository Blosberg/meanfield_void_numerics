//---void.h -function definitions for the void numerics calculation.
// ---last updated on  Tue Nov 26 19:18:48 CET 2013  by  Brendan.Osberg  at location  th-ws-e537

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
t0  = rates_times[0];
t1  = rates_times[1];
rm  = rates_times[2];
double muN_original  = rates_times[3];
E0  = rates_times[4];
CGF = rates_times[5];	//----coarse-graining factor.

L        = sizes[0];
kHNG_CG  = sizes[1];
w        = sizes[2];

shouldplotvoiddist   = B[0];
shouldplotrhos       = B[1];

HNG                  = B[2];
SNG                  = B[3];
LNG                  = B[4];
boltzmann_on_add     = B[5];
boltzmann_on_removal = B[6];
boltzmann_on_uphill  = B[7];

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
attempt=0;

mean 		 = 0.0;
std_dev 	 = 0.0;

//--------------------ALLOCATE AND ASSIGN THE 2-BODY INTERACTION -------------------
xcoarse  = new double[L+1];

//--------------coarse graining here:--------------------------------------

muN_CG = muN_original + gsl_sf_log(CGF); //---scale the effective chemical potential according to the coarse-graining
rp = gsl_sf_exp(muN_CG); 
r  = rp/rm;		//---keep these two lines constant for now!

kHNG_original = 147;

if(SNG || LNG )
	{
	//-----this is still valid regardless of when we TAKE the boltzmann factor.
	//-----we're just determining what it is exactly.

	Bzman_v2 = new double[L+1];
	v2       = new double[L+1];
	p=2*w+1;
	v2_uncoarsened    = new double[p];


	for(i=0;i<=L;i++)
		{
		v2[i]       = 0.0;
		Bzman_v2[i] = 1.0;
		xcoarse[i]  = 0.0;
		}

	for(i=0;i<p;i++)
		{v2_uncoarsened[i]=0.0;}
	if(SNG)
		{
		VNN_SNG_calc( v2_uncoarsened, w, E0);
		}
	else if(LNG)
		{
		VNN_LNG_calc( v2_uncoarsened, w, E0);
		}


	coarse_grain( v2_uncoarsened, xcoarse, v2, L, 2*w+1, CGF); //--coarse-grain the 2-bo	

	//-------------------- got the potential and coarse-grained it: from here, SNG and LNG are treated the same 
	//-------------------- (they just have a different v2 potential from here).

	for(i=0;i<=L;i++)
		{
		Bzman_v2[i] = gsl_sf_exp(-v2[i]);
		}
	}
else if (HNG)
	{
	Bzman_v2 = NULL;	//----- code it this way to act as a warning in case v2 
	v2       = NULL;	//----- is ever accessed during an HNG run (it shouldn't be).
	v2_uncoarsened    = NULL;

	/*-------here is where we hard code the HNG in the same way as SNG.------
	for(i=0;i<=L;i++)
		{
		if(i<kHNG_CG)
			{
			Bzman_v2[i] = 1.0;
			}		
		else if(i>=kHNG_CG)
			{
			Bzman_v2[i] = 0.0;
			}		
		}
	-------------------------------------------------------------------------*/
	}



//-------------------- logarithmic spacing of observation points -------------------------

total_obs_vdist = ceil(   Np10* log10( t1/t0)  );

printq        = new bool[total_obs_vdist];
tpoints_vdist = new double[total_obs_vdist];

for(i=0;i<total_obs_vdist;i++)
	{
 	printq[i]=false;
	tpoints_vdist[i]=0.0;
	}

space_tpoints_logarithmic(t0, t1, Np10 , total_obs_vdist, tpoints_vdist );
//----------------------------------------------------------------------------------------


rho     = 0.0;
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
	delete [] v2_uncoarsened;
	}

delete [] printq;
delete [] tpoints_vdist;

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
//***************************************************************************
bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;    
        (void) closedir (pDir);
    }

    return bExists;
}
//****************************************************************
ODEdat& ODEdat::operator= (const ODEdat& param)
{
v2          	= param.v2;		//----- v2[x]       = the energetic potential of two particles with 'x' empty sites between them.
Bzman_v2	= param.Bzman_v2;	//----- Bzman_v2[x] = exp(-v2[x]), the boltzmann factor thereof.
size_v2		= param.size_v2;		// the number of points in the v2 potential (beyond which it is assumed v2=0)

rho		= param.rho;		//---the average density.
rhostar		= param.rhostar;	//---the normalization factor.
//! rhodot		= param.rhodot;		//---the derivative
empty_space	= param.empty_space;
C1		= param.C1;
C2		= param.C2;

rm		=param.rm; // = r-'minus' -the off rate.
rp		=param.rp; // = r-'plus'  -the on rate.
r		=param.r;  // = r+/r-

kHNG_CG		=param.kHNG_CG; //---the finite size of the particles.
w		=param.w; //---the range of interaction.
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
	mean += i*CGF*(V[i]/Voidsum); //i*CGF is the void size V[i]/Voidsum is the probability.
	checknorm += V[i]/Voidsum; //---normalization factor for Void distances;
	}

for(i=0;i<=L;i++)
	{
	std_dev += ( (V[i]/Voidsum)* (double(i)*CGF-mean)*(double(i)*CGF-mean) ); //---calculate the std.dev.
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

string detail1;
string detail2;

if(HNG)
	{
	detail1="HNG_k";
	detail2=convertInt(kHNG_CG);
	}
else
	{
	detail1="SNG_w";
	detail2=convertInt(w);
	}
//----------------------------------------------------------------
if( plotnum < 10)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_0000%d.txt", pathout.c_str(),plotnum);
	}
else if(plotnum < 100)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_000%d.txt", pathout.c_str(),plotnum);
	}
else if(plotnum < 1000)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_00%d.txt", pathout.c_str(),plotnum);
	}
else if(plotnum < 10000)
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_0%d.txt", pathout.c_str(),plotnum);
	}
else 
	{
	sprintf(cpath, "%sVdist_x_VpL_Vprob_%d.txt", pathout.c_str(),plotnum);
	}

ofstream fvout(cpath);

/*-----------NO NEED TO BOTHER WITH THIS SHIT ANYMORE----------------
//-------get Redner's formula:
double A     = coverage*coverage*(CGF*L) /( kHNG_original * (coverage + kHNG_original*(1.0-coverage )  ) );
double alpha = gsl_sf_log( 1 + coverage/(kHNG_original*(1.0-coverage) )  );
double A_CG     = coverage*coverage*(L) /( kHNG_CG * (coverage + kHNG_CG*(1.0-coverage )  ) );
double alpha_CG = gsl_sf_log( 1 + coverage/(kHNG_CG*(1.0-coverage) )  );
*/

//--------------------------PLOT THE VOID DISTRIBUTION AS A FUNCTION OF X ------------------------
for(i=0;i<L;i++)
	{
	fvout << (i+1)*CGF << " \t " << V[i]/(CGF*L) << " \t " << V[i]/rho << endl;
	}
//---------------NOW THE TIME ASSOCIATED WITH THIS PARTICULAR plotnum iteration------------------

sprintf(cpath, "%splotnum_t_rho_coverage.txt", pathout.c_str());
ofstream pn;

if(plotnum==0)
	pn.open(cpath);
else
	pn.open(cpath,fstream::app);

pn << plotnum << " \t " << t << " \t " << (rho/(L*CGF)) << " \t " << coverage << endl;

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
int  kn = P->kHNG_CG; // ----the particle size.
int  L  = P->L; // ----the system size.
//---------------------------------------------------------
P->attempt +=1;		// keep track of the number of times we've "attempted to time-step here"
P->t = t;		// the current time.

double rm = P->rm;
double rp = P->rp;	// rp/m = 'plus/minus' --the on off rates.

//--------------------  initialize  ----------------------

for(j=0; j<=L ;j++)	
	{ 
	f[j] = 0.0;
	}


(*P).get_rho_anal(V);  //---returns exactly the number of particles in the system -WE DON'T WORRY ABOUT COARSE-GRAINING HERE!
(*P).get_rhostar(V); 
(*P).get_empty_space(V);
 
//----make sure to update rho to fit the numerics in the convolution term (number 4)

//============= here handle finite-size special cases =================================================

//----completely empty system  		0|---------------------------------|L	
f[L]= -1*rp*L*V[L] + rm*V[L-kn] ;	//--only terms 3,4 survive here.


//----single particle in system         0|------------|<-kn->|--------------|L	
f[L-kn]  = -1*rm*V[L-kn] + rp*L*V[L];     //---term 1 only has a single rm-removal rate., 
					// term 2 has full range of 'L' binding sites.

if ( (L-2*kn+1) >=1)
	f[L-kn] += -rp*(L-2*kn+1)*V[L-kn];    //---term 3 : same as canonical case.

//----term (4) ---creation of void 'm' from desorption of neighbour with next void adding up to m.   

for(j=0; j<= (L-2*kn);j++)	 
	{ 
	if( !isnan(1.0/(P->rhostar)) && !isinf(1.0/(P->rhostar)) ) //---- float-error catch (in case rho=0)
		{
		f[L-kn] += (rm/(P->rhostar))*V[j]*V[L-2*kn-j];
		}
	}
//============ done handling special finite-size cases; the rest are canonical ========================
for(m=0; m<= (L-2*kn) ;m++)
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
//---------------------------------------------------------
P->attempt +=1;		// keep track of the number of times we've "attempted to time-step here"
P->t = t;		// the current time.

double rm = P->rm;
double rp = P->rp;	// rp/m = 'plus/minus' --the on off rates.

//--------------------  initialize  ----------------------

for(j=0; j<=L ;j++)	
	{ 
	f[j] = 0.0;
	}

(*P).get_rho_anal(V);  //---returns exactly the number of particles in the system
(*P).get_rhostar(V); 
(*P).get_empty_space(V);
 
//----make sure to update rho to fit the numerics in the convolution term (number 4)

if( P->boltzmann_on_add)
	{
	//======= here handle finite-size special cases =========================

	//----completely empty system  		0|---------------------------------|L	
	f[L]= -1*rp*L*V[L]*P->Bzman_v2[L-1] + rm*V[L-1] ; //--only terms 3,4 survive here.

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

	//---- done handling special finite-size cases; the rest are canonical ==========
	for(m=0; m< (L-1) ;m++)
		{
		//----term(1) -particles on either side leaving.
		f[m]=-2*rm*V[m]; 
		//----term (2) ---creation of void 'm' from a larger one.  
		for(i=m+1; i<=(L-1) ;i++)	
			{ //---term 2 
			f[m] += 2*rp*V[i] * ( (P->Bzman_v2[m]*P->Bzman_v2[i-1-m] )/ (P->Bzman_v2[i]) ); 
			}

		
		//---term(3) disappearance of V[m] due to adsorption in its interior.	
		for(i=0; i<=(m-1) ;i++)	
			{ //---term 3----
			f[m] += -rp*V[m] * ( (P->Bzman_v2[i]*P->Bzman_v2[m-1-i] )/ (P->Bzman_v2[m]) ); 
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

	// N.B. the Bzman array is always exp(-v2(x)), and therefore always <=1.
	// -in order to speed up removal, we have to DIVIDE rates by it.

	//---- L ------------------------------------------	
	f[L]= -1*rp*L*V[L] + rm*V[L-1]/(P->Bzman_v2[L-1]) ; //--only terms 3,4 survive here.

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

	//---- done handling special finite-size cases; the rest are canonical ==========
	for(m=0; m< (L-1) ;m++)
		{
		//----term(1) -particles on either side leaving.
		for(i=0; i<=(L-2-m) ;i++)	
			{
			f[m] -= 2*rm*V[m] * (V[i]/P->rhostar) * (1/(P->Bzman_v2[m]*P->Bzman_v2[i] ))*P->Bzman_v2[i+m+1] ; 
			}
		//----term (2) ---creation of void 'm' from a larger one.  

		for(i=m+1; i<=(L-1) ;i++)	
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
	for (n=0;n<=(L-3*kHNG_CG);n++)
		{
		inner_sum = 0.0;

		for (j=0;j<=n;j++)
			{
			inner_sum += (V[j]*V[n-j]);
			}

		num += (double(n)+2.0*double(kHNG_CG))*inner_sum;
		}
	n=L-2*kHNG_CG;//=1
	inner_sum = 0.0;

	for (j=0;j<=n;j++)
		{
		inner_sum += (V[j]*V[n-j]);
		}

	num += (double(n)+2.0*double(kHNG_CG))*inner_sum;

	//------------GET THE HNG DENOMINATOR ----------
	for (m=0;m<=(L-2*kHNG_CG);m++)
		{
		denom +=2*(kHNG_CG+m)*V[m];
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
	coverage = rho*kHNG_CG/L;
	}
else if(SNG || LNG )
	{
	coverage = rho*p/L; //--implement the "coverage from DNA perspective later.
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
for(m=0;m<=(L-2*kHNG_CG);m++)
	{
	C1 += (-2*rm)*(kHNG_CG+m)*V[m];
	}
//----------------------------------------------
for(n=0;n<=(L-3*kHNG_CG);n++)
	{
	inner_sum=0.0;
	
	for(j=0;j<=n;j++)
		{
		inner_sum +=V[j]*V[n-j];		
		}
	C2+= (n+2*kHNG_CG)*rm*inner_sum;
	}
//---------now add on the last one for L-2kHNG_CG-----------
n=L-2*kHNG_CG;
inner_sum = 0.0;

for (j=0;j<=n;j++)
	{
	inner_sum += (V[j]*V[n-j]);
	}
C2+= (n+2*kHNG_CG)*rm*inner_sum;
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