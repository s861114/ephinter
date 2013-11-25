#ifndef WKSP_H
#define WKSP_H

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <omp.h>
#include <gsl/gsl_sf_trig.h>

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
	double hv_a;
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

	int N_q;
	double qc;
	double h_q;

	int N_phiq;
	double h_phiq;

	int N_omg;
	double omgc;
	double h_omg;
	//----------------------------------



	
	
	//-----------Variables---------------
	int vg;
	int numofthread;
	int N;  // num of layer
	int N2; // dim of system.
	int tot_th;
	// N2= 2*N for full band model, N2 = 2 for 2band model

	int q_idx;
	int phiq_idx;
	int omg_idx;

	int k_idx;
	int phik_idx;

	double**** polar_total_var; 
	double****** polar_band_var;
	double**** polar_total_var_analytic; 

	double q_real;
	double phiq_real;
	double omg_real;

	double Ef;
	//--------------------------------------






	//-----functions&associated varables---
	void initial_define_constant(void);

	void initial_stacking_determiner(void);
	
	void band_cal(void);

	void band_cal_q(void);
	
	void polar_func(void);

	void sub_polar_func(void);

	void sum_total_func(void);

	void sum_total_func_analytic(void);

	double overlap(int i, int j);

	void print_polar(void);

	double Re_polar_analytic_1(double x, double v);

	double Re_polar_analytic_2(double x, double v);

	double f1(double x, double y);
	
	double f2(double x, double y);

	double f3(double x, double y);

	double signum(double x);
	//-------------------------------------





	/// memory
	void initial_malloc(void);
	gsl_eigen_hermv_workspace** ws;
	gsl_vector** eval;

	gsl_matrix_complex*** eigen_state; //matrix* eigenstate[radial][theta]
	gsl_matrix_complex*** eigen_state_q; //matrix* eigenstate[radial][theta]
	double*** energy; // energy[band][radial][theta]
	double*** energy_q; // energy[band][radial][theta]
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
