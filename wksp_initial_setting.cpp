#include "wksp.h"
#include "function.h"

WKSP::WKSP()
{
	initial_define_constant();
	initial_read_setting();
	set_H();
}

void WKSP::initial_define_constant(void)
{

	gamma0=3; // 3eV
	gamma1=0.3;//0.3eV
	a=2.46; // lattice constant 2.46 angstrom
	d=3.35; // distance between layers 3.35 angstrom
	hv_a=sqrt(3)/2*gamma0*a;
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
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&kc);		
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_theta);		
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_idf);	

	h_radial = kc / N_radial;
	h_theta = 2.0*M_PI/ N_theta;
	c_theta = 1.0/N_theta;	

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
	H = mat_gsl_matrix_complex_pointer.Tmatrix2(N_radial,N_theta);
	epho = mat_gsl_matrix_complex_pointer.Tmatrix3(N2,N_radial,N_theta);
	// N2 is # of bnads of system.
	energy = mat_double.Tmatrix3(N2,N_radial,N_theta);
	en = mat_double.Tmatrix1(N2);

	for(int i=0; i<N_radial; i++)
	{
		for(int j=0; j<N_theta; j++)
		{
			eigen_state[i][j]= gsl_matrix_complex_alloc(N2,N2);
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
	mat_gsl_matrix_complex_pointer.Tfree2(H);
	mat_double.Tfree3(energy);
	mat_double.Tfree1(en);
	printf("finish memory free\n");
}


