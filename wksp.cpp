#include "wksp.h"
#include "function.h"
#include <cmath>
void WKSP::band_cal(void)
{
	FILE *fp=fopen("/media/sf_LinuxShare/output.txt","w");
	omp_set_num_threads (3);
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		
		int cur_th = omp_get_thread_num();
		fprintf(stderr,"cur_th : %d\n", omp_get_num_threads());
		// getting current thread num
		for(int j=0; j<N_theta; j++)
		{
			gsl_eigen_hermv(H[i][j],eval[cur_th],eigen_state[i][j],ws[cur_th]);
			gsl_eigen_hermv_sort(eval[cur_th],eigen_state[i][j],GSL_EIGEN_SORT_VAL_ASC);

			for(int e=0; e<N2; e++)
				energy[e][i][j] = eval[cur_th]->data[e];			
		}
	}
				
	for(int i=0;i<N_radial;i++)
	{
		fprintf(fp,"%f ",i*h_radial);
		for(int e=0;e<N2;e++)
		{	
			fprintf(fp,"%f ",energy[e][i][0]);
		}
		fprintf(fp,"\n");
	}


}

void WKSP::band_cal_q(void)
{

//	FILE *fp=fopen("/media/sf_LinuxShare/output.txt","w");
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		int cur_th = omp_get_thread_num();
		// getting current thread num
		for(int j=0; j<N_theta; j++)
		{
			gsl_eigen_hermv(H[i][j],eval[cur_th],eigen_state_q[i][j],ws[cur_th]);
			gsl_eigen_hermv_sort(eval[cur_th],eigen_state_q[i][j],GSL_EIGEN_SORT_VAL_ASC);
			
			for(int e=0; e<N2; e++)
				energy_q[e][i][j] = eval[cur_th]->data[e];			
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
	for (q_idx=1; q_idx<N_q; q_idx++)
	{
		for (phiq_idx=0; phiq_idx<N_phiq; phiq_idx++)
		{
			for (omg_idx=0; omg_idx<N_omg; omg_idx++)
			{
				fprintf(stderr,"%d  /  %d\n",q_idx,N_q);
				q_real=q_idx*h_q;
				phiq_real=phiq_idx*h_phiq;
				omg_real=omg_idx*h_omg;
				set_H();
				band_cal_q();				
				sub_polar_func();
			}
		}
	}
	sum_total_func(); //This function sums sub_polar_var into polar_total_var
}



void WKSP::initial_stacking_determiner(void)
{
	switch(N_idf)
	{
		case 1:N=1;N2=2;break;
		case 2:N=2;N2=4;break;
	}
}


void WKSP::sub_polar_func(void)
{
	//Real Part
	const double constant=-1/M_PI/M_PI/hv_a/a;
	double k_real;
	double en_copy;
	double en_q_copy;
	double sum=0;
	double denom;
	for (int bnd_idx1=0; bnd_idx1<N2; bnd_idx1++)
	{
		for (int bnd_idx2=0; bnd_idx2<N2; bnd_idx2++)
		{

			//-----------------Real Integration Part (CORE)---------------------
			sum=0;			
			for (k_idx=0;k_idx<N_radial; k_idx++)
			{
				for (phik_idx=0;phik_idx<N_theta;phik_idx++)
				{
					k_real		=k_idx*h_radial;

					en_copy		=  energy[bnd_idx1][k_idx][phik_idx];
					en_q_copy	=energy_q[bnd_idx2][k_idx][phik_idx];
					if ( (Ef>en_copy)&&(Ef<en_q_copy)==true ||  (Ef<en_copy)&&(Ef>en_q_copy)==true )
					{
						denom=(omg_real+en_copy-en_q_copy);
						if(denom!=0)
							sum+=overlap(bnd_idx1, bnd_idx2)*h_radial*k_real/denom; //Here integration is performed.					
						if ( (Ef>en_copy)&&(Ef<en_q_copy) ) 
						   sum=-sum;	
					}
				}
			}
			polar_band_var[q_idx][phiq_idx][omg_idx][bnd_idx1][bnd_idx2][0]=sum*constant;

			//-------------------------------------------------------------------

			//-----------------Imaginary Integration Part (CORE)-----------------
			sum=0;	
			for (int k_idx=0;k_idx<N_radial; k_idx++)
			{
				for (int phik_idx=0;phik_idx<N_theta;phik_idx++)
				{
//					sum+=1; //Here integration is performed.
				}
			}
			const double constant=-1/M_PI/M_PI/hv_a/a;
			polar_band_var[q_idx][phiq_idx][omg_idx][bnd_idx1][bnd_idx2][0]=sum*constant;
			//-------------------------------------------------------------------
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


double WKSP::overlap(int i, int j)
{
	for(int bnd_idx=0;bnd_idx<N2;bnd_idx++)
	{
		gsl_complex temp1=gsl_complex_conjugate(	gsl_matrix_complex_get(eigen_state[k_idx][phik_idx],bnd_idx,i)	);
		gsl_complex temp2=	gsl_matrix_complex_get(eigen_state_q[q_idx][phiq_idx],bnd_idx,j);
		return	GSL_REAL(gsl_complex_mul(temp1,temp2));
	}	
}
