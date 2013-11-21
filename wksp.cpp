#include "wksp.h"

void WKSP::band_cal(void)
{
	FILE *fp=fopen("/media/sf_LinuxShare/output.txt","w");
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		int cur_th = omp_get_thread_num();
		// getting current thread num
		for(int j=0; j<N_theta; j++)
		{
			gsl_eigen_hermv(H[i][j],eval[cur_th],eigen_state[i][j],ws[cur_th]);
			gsl_eigen_hermv_sort(eval[cur_th],eigen_state[i][j],GSL_EIGEN_SORT_VAL_ASC);
			
			for(int e=0; e<N2; e++)
				energy[e][i][j] = eval[cur_th]->data[e];			
		}
	}
/*				
	for(int i=0;i<N_radial;i++)
	{
		fprintf(fp,"%f ",i*h_radial);
		for(int e=0;e<N2;e++)
		{	
			fprintf(fp,"%f ",energy[e][i][0]);
		}
		fprintf(fp,"\n");
	}
*/

}
