//---------------------------------------------------------------------------------------------
// ---last updated on  Sun Apr 27 21:18:49 CEST 2014  by  bren  at location  , bren-Desktop

//  changes from  Sun Apr 27 21:18:49 CEST 2014 : implemented criteria to terminate the run when a maximum in the density is encountered -called "should_peakterminate"

//  changes from  Fri Apr 25 17:08:48 CEST 2014 : optimized the process by iteratively cutting off large voids that have gone negative. cutoff is in func_SLNG"

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
double rho_t=0.0, rho_old=0.0, rho_old2=0.0;
double t_old=0.0, t_old2=0.0, rhodot_anal=0.0, rhodot_num=0.0, rhodot_num_old=0.0, rhodot_num_old2=0.0;
//----------set up /initialize the array --------------
double V[NV];
for(i=0;i<NV;i++)
	V[i]=P->V_IC[i]; //---the constructor sets up the initial V array, and then V_IC is not used.
			 //--- without an import setting, this is just all zeros except for at 

double t = P->t; //--- likewise 't' value may be imported from an unfinished run. 
		 //--- Either way it's set in the P constructor.

//---------------------------------------------------------------------------- 

//--------------  document the potential in the log file.  ------------------------
*P->log << "\n-------------------------------------------"; 

bool shouldprint=false;
//--------------SET UP THE NECESSARY dE evolution parameters--------------	
const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;

double h = P->t0;	// t0 is not changed from the input file in the constructor 
			// regardless of IC condition. This is just how large we step initially.
gsl_odeiv_step    * s      =  gsl_odeiv_step_alloc (T, NV);
gsl_odeiv_control * con    =  gsl_odeiv_control_y_new (h,h);
gsl_odeiv_evolve  * e      =  gsl_odeiv_evolve_alloc (NV);


gsl_odeiv_system sys_HNG = {func_HNG, jac, NV, P};
gsl_odeiv_system sys_SLNG = {func_SLNG, jac, NV, P};

int status=0;
//*****************   NOW  TIME STEP **********************

(*P->log) << "\n just before beginning time evolution, t1 is =" << P->t1 << endl;

long int step=0;
double empty_space;

const double Np10_rho =30; //---this is the resolution cut-off for plotting FOR THE DENSITY CURVE
			   //--- IT HAS A DIFFERENT MEANING FROM THE Np10 that we use for the void dists.

const double deltaspacing = gsl_sf_log(10)/Np10_rho; //steps with lower resolution than this don't get plotted


rho_t	   = (*P).get_rho_anal(V);
P->step    = 0;

double t_lastplot = t+h/10.0; //---this should ensure that the first point gets plotted

while (t < P->t1)
	{ 
	P->t = t;

	t_old2   = t_old;
	rho_old2 = rho_old; //---update iteratively.
	
	t_old = t;
	rho_old =  (*P).get_rho_anal(V);//=rho_t;

	rhodot_num_old2 = rhodot_num_old;
	rhodot_num_old  = rhodot_num;
	
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
	
	
	//--------- check which V's have gone negative and eliminate them for efficiency -----------
	
	if ( P->should_check_neg )
		{

		while( V[P->phys_bound] < 0.0 ) 
			{
			P->has_been_neg[P->phys_bound] = true;
			V[P->phys_bound] 	       =  0.0; //--- hard set it to zero and don't worry about it any more. the f's will all be zero too.
			P->phys_bound--;
	
			//--------only relevant for HNG case -----------
			if(P->HNG)
				{
				if( P->phys_bound > (P->L-P->a) && P->phys_bound < (P->L) )
					{
					P->phys_bound = (P->L-P->a);
					for(i=(P->L-P->a+1); i<(P->L); i++)
						{
						P->has_been_neg[i] = true;
						}
					}

				if( P->phys_bound > (P->L)-2*(P->a) && P->phys_bound < ((P->L)-(P->a)) )
					{
					P->phys_bound = (P->L)-(2*(P->a));
					for(i=( (P->L)-2*(P->a)+1); i<(P->L)-(P->a); i++)
						{
						P->has_been_neg[i] = true;
						}
					}

				}
			//-----------------    DOWN TO HERE   -----------

			/*-------
			if(P->phys_bound < 2*P->a)
				{
				cout << "\n ERROR: phys_bound = " << P->phys_bound << " is getting suspiciously small at t=" << P->t << endl;
				exit(1);
				}
			-------*/
			}
		}

	//------------------------- DONE ELIMINATING GARBAGE NEGATIVE V'S.  NOW PLOT rho v time ----------------------

	if( P->shouldplotrhos && ( gsl_sf_log(t/t_lastplot) > deltaspacing ) ) 
		{
		// output to the t vs. rho file.
		t_lastplot = t; //---update the current "last point plotted"
		rho_t	    = (*P).get_rho_anal(V);
		dummy       = (*P).get_mean_stddev(V);
		//---------------------OUTPUT TO FILE ----THIS IS WHERE WE PLOT THE FILLING -------------------------
		 *foutmain << t << "\t" << ((*P).rho)/(P->L ) << "\t" << (rhodot_num)/(P->L ) << " \t " << P->phys_bound << endl;
		 //removed: << "\t" << P->mean << " \t " << P->std_dev << "\t" << rhodot_num << endl;
		 //----- we're not using the mean or std. dev.'s anymore, so don't bother with them. 
		}

	
    if( P->should_peakterminate &&  t > P->t_export) //---dont even think about exiting until after t_export.
	{
    	if( ( (P->rho) <= (rho_old) ) &&  ( (rho_old) <= (rho_old2) )   )
    	{//----two successive steps down or nowhere in density after t_export ==> stopping.

    		*foutmain << t_old2 << "\t" << (rho_old2)/(P->L ) << "\t" << (rhodot_num_old2)/(P->L ) << endl;
    		*foutmain << t_old  << "\t" << (rho_old )/(P->L ) << "\t" << (rhodot_num_old )/(P->L ) << endl;
    		*foutmain << t      << "\t" << (rho_t)/(P->L )    << "\t" << (rhodot_num     )/(P->L ) << endl;
   		 
    		if(P->should_export_IC)
    		{
    			P->export_IC(V);
    		}
    		
    		(*P->log) << "\n break condition encountered. slope is negative.";
	   		(*P->log) << " Exporting current V-distribution, and terminating run at t= " << P->t << endl;
	
	   		P->rho = -777.7*P->L;	//---to make sure we don't accidentally take this value to be an equilibrated termination point.
    		break;
	
    	}
	}
    
    }//================================================================-----------------------------
//finished while loop (i.e. t has reached t1)

if(P->should_export_IC)
	{
	P->export_IC(V);
	}

//--we don't bother with this anymore, for efficiency. (*P->log) << "\n at the end of the simulation, " << step << " steps were taken." << endl;
// cout      << "\n at the end of the simulation, " << step << " steps were taken." << endl;

       gsl_odeiv_evolve_free  (e);
       gsl_odeiv_control_free (con);
       gsl_odeiv_step_free    (s);


return 1;
}


//****************************************************************************

