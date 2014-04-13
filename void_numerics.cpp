//---------------------------------------------------------------------------------------------
// ---last updated on  Sun Apr 13 23:34:52 CEST 2014  by  bren  at location  bren-Desktop

//  changes from  Sun Apr 13 23:34:52 CEST 2014 : modified the script for explicitly-sized particles (i.e. suitable for dimers, trimers, etc., no longer just coarse-grained full-sized Nucl's.)
//---------------------------------------------------------------------------------------------


//---------- IN THIS FILE WE MERELY CALL THE RESPECTIVE func_HNG/SLNG functions 
//---------- and they make the distinction of which way to implement the boltzmann weights.


#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_log.h>



//------------------------------------------
#include "void_numerics.h"
//  #include "void.h"
//  #include <bren_lib.h>
//------------------------------------------



using namespace std;

// int  charlength=400; //--the number of characters in the string for our path.
//************************************************************************************
int jac (double t, const double c[], double *dfdy, double dfdt[], void *params)
{
//-----we only have a single vector here; no F(x1,x2...xN); therefore no Jacobian
}



//***************************************************************

int time_go(ODEdat *P, ofstream * foutmain)
{
int   charlength=400; //--the number of characters in the string for our path.
int   NV = P->L+1; //---NV is the number of voids to be considered (including 0)

char cpath[charlength];
int i,dummy;
double rho_t=0.0, rho_old=0.0, t_old=0.0, rhodot_anal=0.0, rhodot_num=0.0;
//----------set up /initialize the array --------------
double V[NV];
for(i=0;i<=NV;i++)
	V[i]=0.0;

//---------------------------INITIALIZE THE ARRAY ------------------------------
// (*P).initialize_V(init,V); //---init is the size of the voids you want to start with.

V[P->L]=1.0; //----just set the whole system to completely empty for the beginning.

//---------------------------------------------------------------------------- 

//--------------  document the potential in the log file.  ------------------------
*P->log << "\n-------------------------------------------"; 
/*
*P->log << "\n the two-body is of size " << P->size_v2 << ", and is as follows: \n ";

for (i=0;i< P->size_v2;i++)
	{
	*P->log << i << " \t " <<  P->v2[i] << " \t " << P->Bzman_v2[i] << endl;
 	}
*P->log << "\n-------------------------------------------"; 
*/

bool shouldprint=false;
//--------------SET UP THE NECESSARY dE evolution parameters--------------	
const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;

double h = P->t0;
gsl_odeiv_step    * s      =  gsl_odeiv_step_alloc (T, NV);
gsl_odeiv_control * con    =  gsl_odeiv_control_y_new (h,h);
gsl_odeiv_evolve  * e      =  gsl_odeiv_evolve_alloc (NV);


gsl_odeiv_system sys_HNG = {func_HNG, jac, NV, P};
gsl_odeiv_system sys_SLNG = {func_SLNG, jac, NV, P};

int status=0;
//*****************   NOW  TIME STEP **********************

(*P->log) << "\n just before beginning time evolution, t1 is =" << P->t1 << endl;

long int step=0;
double t=0.0;
double empty_space;

const double Np10_rho =30; //---this is the resolution cut-off for plotting
const double deltaspacing = gsl_sf_log(10)/Np10_rho; //steps with lower resolution than this don't get plotted


rho_t	   = (*P).get_rho_anal(V);
P->step    = 0;

double t_lastplot=h/10.0; //---this should ensure that the first point gets plotted

while (t < P->t1)
	{ 
	step ++;
	P->t = t;

	P->attempt = 0;

	t_old = t;
	rho_old =  (*P).get_rho_anal(V);//=rho_t;

	if(P->HNG)
		{
		status = gsl_odeiv_evolve_apply (e, con, s, &sys_HNG, &t, P->t1, &h, V);
		}
	else if( P->SNG || P->LNG )
		{
		status = gsl_odeiv_evolve_apply (e, con, s, &sys_SLNG, &t, P->t1, &h, V);
		}
	else
		{
		cout << "\n NG type undefined.\n";
		}

	if (status != GSL_SUCCESS)
		{               
		cout << "\n process interrupted at t=" << t << ", status failed in Re inc." << endl;
		break;
		}
	P->step ++;
	//------  get analytic and numeric rho, rhodot for comparison-------------

	rho_t	    = (*P).get_rho_anal(V);
	empty_space = (*P).get_empty_space(V);
	dummy       = (*P).get_mean_stddev(V);

	rhodot_num  = (rho_t-rho_old)/(t-t_old);
	P->t = t;

	if (rho_t > P->maxrho)
		{
		P->maxrho = rho_t;	// remember that the *ACTUAL* density still has to be divided by L
		}

	//--------take snapshots in time------------------------------

	shouldprint = (*P).printtime();

	if( (*P).shouldplotvoiddist && shouldprint ) 
		{
		(*P).printout_voiddist(V);
		P->plotnum+=1;
		}
	//----these don't get printed with the same frequency as the filling rates below
	//-----track the probabilities------------

	if( P->shouldplotrhos && ( gsl_sf_log(t/t_lastplot) > deltaspacing ) ) // this is for hig res plots -use this one for the single timeline plot option.
		{
		t_lastplot = t; //---update the current "last point plotted"

		//---------------------OUTPUT TO FILE ----THIS IS WHERE WE PLOT THE FILLING -------------------------
		 *foutmain << t << "\t" << ((*P).rho)/(P->L ) << "\t" << P->mean << " \t " << P->std_dev << "\t" << rhodot_num << endl;

		}


    }//================================================================-----------------------------
//finished while loop (i.e. t has reached t1)

(*P->log) << "\n at the end of the simulation, " << step << " steps were taken." << endl;
cout      << "\n at the end of the simulation, " << step << " steps were taken." << endl;

       gsl_odeiv_evolve_free  (e);
       gsl_odeiv_control_free (con);
       gsl_odeiv_step_free    (s);


return 1;
}


//****************************************************************************

