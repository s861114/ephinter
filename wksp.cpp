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


void WKSP::polar_func(void)
{
	for (q_idx=0; q_idx<N_q; q_idx++)
	{
		for (phiq_idx=0; phiq_idx<N_phiq; phiq_idx++)
		{
			for (omg_idx=0; omg_idx<N_omg; omg_idx++)
			{
				q_real=q_idx*h_q;
				phiq_real=phiq_idx*h_phiq;
				omg_real=omg_idx*h_omg;
				sub_polar_func();
			}
		}
	}
	sum_total_func(); //This function sums sub_polar_var into polar_total_var
}






void WKSP::sub_polar_func(void)
{
	
	//Real Part
	double sum=0;
	for (int bnd_idx1=0; bnd_idx1<N2; bnd_idx1++)
	{
		for (int bnd_idx2=0; bnd_idx2<N2; bnd_idx2++)
		{

			for (int k_idx=0;k_idx<N_radial; k_idx++)
			{
				for (int phik_idx=0;phik_idx<N_theta;phik_idx++)
				{
					sum+=sum; //Here integration is performed.
				}
			}
			polar_band_var[q_idx][phiq_idx][omg_idx][bnd_idx1][bnd_idx2][0]=sum;

			//Imaginary Part
			sum=0;	
			for (int k_idx=0;k_idx<N_radial; k_idx++)
			{
				for (int phik_idx=0;phik_idx<N_theta;phik_idx++)
				{
					sum+=sum; //Here integration is performed.
				}
			}
			polar_band_var[q_idx][phiq_idx][omg_idx][bnd_idx1][bnd_idx2][1]=sum;
		}
	}
}






void WKSP::sum_total_func(void)
{
	for (int idf=0;idf<2;idf++)
	{
		for (int q_idx=0; q_idx<N_q; q_idx++)
		{
			for (int phiq_idx=0; phiq_idx<N_phiq; phiq_idx++)
			{
				for (int omg_idx=0; omg_idx<N_omg; omg_idx++)
				{
					double sum=0;
					for (int bnd_idx1=0; bnd_idx1<N2; bnd_idx1++)
					{
						for (int bnd_idx2=0; bnd_idx2<N2; bnd_idx2++)
						{
							sum=sum+polar_band_var[q_idx][phiq_idx][omg_idx][bnd_idx1][bnd_idx2][idf];
						}
					}
					polar_total_var[q_idx][phiq_idx][omg_idx][idf]=sum;
				}
			}
		}
	}
}
