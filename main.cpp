#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"
#include <time.h>
int main(void)
{
	WKSP wksp;
	wksp.band_cal();
//	wksp.polar_func();
	wksp.sum_total_func_analytic();	
	wksp.print_polar();
	
	return 0;
}
