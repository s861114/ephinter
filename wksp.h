#ifndef WKSP_H
#define WKSP_H

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <omp.h>


class WKSP  // wksp = workspace
{
public:
	WKSP();  // initailize, memory allocation, define constant;
	~WKSP(); // memory free





	//------- constant_setting-----------

	double gamma0;
	double gamma1;
	double a;
	double d;
	//------------------------------------





	//----------Setup.set-----------------
	void initial_read_setting(void);
	
	int N_radial; 
	double kc;
	double h_radial; // kc/N_radial;
	
	int N_theta;
	double h_theta; // 2pi/N_theta;
	double c_theta; // 1/N_theta;

	int N_idf;	
	//----------Setup.set End--------------



	
	
	//-----------Variables---------------
	int vg;
	int numofthread;
	int hv_a;
	int N;  // num of layer
	int N2; // dim of system.
	// N2= 2*N for full band model, N2 = 2 for 2band model
	//--------------------------------------

	//--------------Functions--------------
	void initial_define_constant(void);

	//-------------------------------------

	int tot_th;

	/// memory
	void initial_malloc(void);
	gsl_eigen_hermv_workspace** ws;
	gsl_vector** eval;

	gsl_matrix_complex*** eigen_state; //matrix* eigenstate[radial][theta]
	gsl_matrix_complex**** epho; //density_matrix //matrix* epho[band][radial][theta]
	double*** energy; // energy[band][radial][theta]
	void set_en(void); // calculate electron density
	double* en; // electron density[layer]

	gsl_matrix_complex*** H; //matrix* Hamailtonian[radial][theta]

//------------------Build Hamiltonians--------------------
	void set_H(void);

	void set_H_A(void);
	//2band model
	void set_H_AB2(void);
	void set_H_ABC2(void);
	void set_H_ABCA2(void);

	// full band model
	void set_H_AB(void);
	void set_H_ABC(void);
	void set_H_ABCA(void);
	void set_H_ABCAB(void);

	void set_H_ABA(void);
	void set_H_ABAB(void);
	void set_H_ABABA(void);
//--------------------------------------------------


};

#endif
