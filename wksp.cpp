#include "wksp.h"
#include "function.h"
#include <cmath>
void WKSP::band_cal(void)
{
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
		fprintf(stderr,"%d  /  %d\n",q_idx,N_q);
		for (phiq_idx=0; phiq_idx<N_phiq; phiq_idx++)
		{
			for (omg_idx=0; omg_idx<N_omg; omg_idx++)
			{
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
			polar_band_var[q_idx][phiq_idx][omg_idx][bnd_idx1][bnd_idx2][1]=sum*constant;
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

void WKSP::sum_total_func_analytic(void)
{
	double conv_q=0;
	double conv_w=0;
	for (int q_idx=0; q_idx<N_q; q_idx++)
	{
		for (int phiq_idx=0; phiq_idx<N_phiq; phiq_idx++)
		{
			for (int omg_idx=0; omg_idx<N_omg; omg_idx++)
			{
				q_real = q_idx * h_q;
				omg_real = omg_idx*h_omg;
				conv_q = q_real * hv_a / (a*Ef);
				conv_w = hv_a / Ef / a * omg_real;
				if (conv_w>conv_q)
					polar_total_var_analytic[q_idx][phiq_idx][omg_idx][0]=Re_polar_analytic_1(conv_q,conv_w);
				else if (conv_w<conv_q)
					polar_total_var_analytic[q_idx][phiq_idx][omg_idx][0]=Re_polar_analytic_2(conv_q,conv_w);
				else
					polar_total_var_analytic[q_idx][phiq_idx][omg_idx][0]=0;
				
			}
		}
	}
}


double WKSP::signum(double x)
{
	if (x>0)
		return 1;
	else
		return -1;
}

double WKSP::Re_polar_analytic_1(double x, double v)
{
	fprintf(stderr,"%f %f, %f\n",x,-v,f2(x,-v));
	double temp;
	temp=1-1/( 8*sqrt(v*v-x*x) )*( f1(x,v)*(fabs(2+v)>x)+signum(v-2+x)*f1(x,-v)*(fabs(2-v)>x)+f2(x,v)*( (x+2>v)+((2-x-v)>0)));
	return temp;
}

double WKSP::Re_polar_analytic_2(double x, double v)
{
	double temp;
	temp=1-1/( 8*sqrt(v*v-x*x) )*( f3(x,v)*(fabs(2+v)<x)+f3(x,-v)*(fabs(2-v)<x) + M_PI*x*x/2*( (abs(v+2)>x) + (abs(v-2)>x) ) );
	return temp;
}

double WKSP::f1(double x, double v)
{
	double temp;
	temp=(2+v)*sqrt( (2+v)*(2+v) - x*x )- x*x * log ( ( sqrt((2+v)*(2+v)-x*x) + (2+v) ) / ( abs( sqrt(v*v-x*x) + v) ) );
	return temp;
}

double WKSP::f2(double x, double v)
{
	double temp;
	temp=x*x*log( (v-sqrt(v*v-x*x))/x );
	return temp;
}

double WKSP::f3(double x, double v)
{
	double temp;
	temp=(2+v)*sqrt(x*x-(2+v)*(2+v))+x*x*asin((2+v)/x);
	return temp;
}


double WKSP::overlap(int i, int j)
{
	double sum=0;
	for(int bnd_idx=0;bnd_idx<N2;bnd_idx++)
	{
		gsl_complex temp1=gsl_complex_conjugate(	gsl_matrix_complex_get(eigen_state[k_idx][phik_idx],bnd_idx,i)	);
		gsl_complex temp2=	gsl_matrix_complex_get(eigen_state_q[q_idx][phiq_idx],bnd_idx,j);
		sum+=GSL_REAL(gsl_complex_mul(temp1,temp2));
	}	
	return sum;
}


void WKSP::print_polar(void)
{
	FILE *fp_num=fopen("/media/sf_LinuxShare/output_num.txt","w");
	FILE *fp=fopen("/media/sf_LinuxShare/output.txt","w");
	for (int idf=0;idf<2;idf++)
	{
		for (int phiq_idx=0; phiq_idx<N_phiq; phiq_idx++)
		{
			for (int omg_idx=0; omg_idx<N_omg; omg_idx++)
			{
				for (int q_idx=0; q_idx<N_q; q_idx++)
				{
					fprintf(fp_num,"%f	%d	%d	%d	%d\n",polar_total_var[q_idx][phiq_idx][omg_idx][idf],q_idx,phiq_idx,omg_idx,idf);
					fprintf(fp,"%f	%d	%d	%d	%d\n",polar_total_var_analytic[q_idx][phiq_idx][omg_idx][idf],q_idx,phiq_idx,omg_idx,idf);
				}
			}
		}
	}
}
