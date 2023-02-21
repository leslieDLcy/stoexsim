#include <cstdio>
#include <ctime>
#include <cstdlib>
#include "RAN3.C"

int main(int argc, char *argv[])
{
	float ndev=0., min=1., max=-1.;
	long int idum;
	if(argc<3)
	{
		fprintf(stderr, "\n\nusage: ndev <mean> <stdev> $RANDOM [min] [max]\n\n");
		exit(0);
	}
	sscanf(argv[3],"%d", &idum);
	if(argc>5){
                sscanf(argv[4],"%f", &min);
                sscanf(argv[5],"%f", &max);
		if(min>max){
			fprintf(stdout, "Minimum is greater than maximum!\n");
			fprintf(stderr, "\n\nusage: ndev <mean> <stdev> [seed] [min] [max]\n\n");
			exit(0);
		}
        }

	if(idum>0)
		idum*=-1;

	ndev = (2.*atof(argv[2])*(ran3(&idum)-0.5))+atof(argv[1]);
	if(min<max){
		while ( ndev < min || ndev > max )
		{
			ndev = (2.*atof(argv[2])*(ran3(&idum)-0.5))+atof(argv[1]);
		}	
	}
	fprintf(stdout, "%f\n", ndev);
	return(0);
}
