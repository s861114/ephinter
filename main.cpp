#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"

int main(void)
{
	WKSP wksp;

	wksp.set_H_AB();
	wksp.calcul_pho();

	return 0;
}
