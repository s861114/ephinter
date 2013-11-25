#include "wksp.h"
#include "function.h"

WKSP::WKSP()
{
	initial_define_constant();
	initial_read_setting();
	initial_stacking_determiner();
	initial_malloc();
	q_real=0;
	phiq_real=0;
	set_H();
}

void WKSP::initial_define_constant(void)
{

	gamma0=3; // 3eV
	gamma1=0.3;//0.3eV
	a=2.46; // lattice constant 2.46 angstrom
	d=3.35; // distance between layers 3.35 angstrom
	hv_a=sqrt(3)/2*gamma0*a; //hbar*velocity
	vg=0;

	// tot_th means num of logical CPU.
	//  = physical CPU # * 2 (hyperthreading)
	#pragma omp parallel
	{
		tot_th = omp_get_num_threads();
	}
}
	
void WKSP::initial_read_setting(void)
{
	FILE* fp = fopen("setup.set","r");
	char dump[100];
	char dump1[100];
	/*while(!feof(fp))
	{
		fgets(dump,100,fp);
		printf("%s",dump);
	}*/

	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_radial);
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&kc);		// k*a unit
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_theta);		
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_idf);	
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_q);	
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&qc);	
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_phiq);	
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_omg);	
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&omgc); // w*a/vF unit	
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&Ef); // w*a/vF unit

	h_radial = kc / N_radial;
	h_theta = 2.0*M_PI/ N_theta;
	h_q = qc / N_q;
	h_phiq = 2.0*M_PI / N_phiq;
	h_omg = omgc / N_omg;		

}


void WKSP::initial_malloc(void)
{
	printf("start initial_malloc()\n");
	umat<gsl_eigen_hermv_workspace*> mat_gsl_eigen_hermv_workspace_pointer;
	umat<gsl_matrix_complex*> mat_gsl_matrix_complex_pointer;
	umat<gsl_vector*> mat_gsl_vector_pointer;
	umat<double> mat_double;


	ws = mat_gsl_eigen_hermv_workspace_pointer.Tmatrix1(tot_th);
	eval = mat_gsl_vector_pointer.Tmatrix1(tot_th);
	for(int i=0; i<tot_th; i++)
	{
		ws[i] = gsl_eigen_hermv_alloc(N2);
		eval[i] = gsl_vector_alloc(N2);
	}
	eigen_state = mat_gsl_matrix_complex_pointer.Tmatrix2(N_radial,N_theta);
	eigen_state_q = mat_gsl_matrix_complex_pointer.Tmatrix2(N_radial,N_theta);
	H = mat_gsl_matrix_complex_pointer.Tmatrix2(N_radial,N_theta);
	// N2 is the number of bands of the system.
	energy = mat_double.Tmatrix3(N2,N_radial,N_theta);
	energy_q = mat_double.Tmatrix3(N2,N_radial,N_theta);
	en = mat_double.Tmatrix1(N2);
	
	polar_band_var=mat_double.Tmatrix6(N_q,N_phiq,N_omg,N2,N2,2);
	polar_total_var=mat_double.Tmatrix4(N_q,N_phiq,N_omg,2);
	polar_total_var_analytic=mat_double.Tmatrix4(N_q,N_phiq,N_omg,2);

	for(int i=0; i<N_radial; i++)
	{
		for(int j=0; j<N_theta; j++)
		{
			eigen_state[i][j]= gsl_matrix_complex_alloc(N2,N2);
			eigen_state_q[i][j] = gsl_matrix_complex_alloc(N2,N2);
			H[i][j]= gsl_matrix_complex_alloc(N2,N2);
		}
	//	printf("malloc %d\n",i);
	}
	printf("finish initial_malloc()\n");

}


WKSP::~WKSP()
{
	printf("start memory free\n");
	umat<gsl_eigen_hermv_workspace*> mat_gsl_eigen_hermv_workspace_pointer;
	umat<gsl_matrix_complex*> mat_gsl_matrix_complex_pointer;
	umat<gsl_vector*> mat_gsl_vector_pointer;
	umat<double> mat_double;
	

	for(int i=0; i<tot_th; i++)
	{
		gsl_eigen_hermv_free(ws[i]);
		gsl_vector_free(eval[i]);
	}

	for(int i=0; i<N_radial; i++)
	{
		for(int j=0; j<N_theta; j++)
		{
			gsl_matrix_complex_free(eigen_state[i][j]);
			gsl_matrix_complex_free(H[i][j]);
		}
	}
	mat_gsl_matrix_complex_pointer.Tfree2(eigen_state);
	mat_gsl_matrix_complex_pointer.Tfree2(eigen_state_q);
	mat_gsl_matrix_complex_pointer.Tfree2(H);
	mat_double.Tfree3(energy);
	mat_double.Tfree3(energy_q);
	mat_double.Tfree1(en);
	mat_double.Tfree4(polar_total_var);
	mat_double.Tfree6(polar_band_var);
	printf("finish memory free\n");
}


