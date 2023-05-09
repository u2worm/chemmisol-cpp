#include "chemical.h"

int main(int, char *[])
{
	Problem problem;
	problem.distanceToEquilibrium();
	problem.setPH(7.6);
	problem.distanceToEquilibrium();
	
	return 0;
}
