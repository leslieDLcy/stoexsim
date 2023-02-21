#include <cstdio>
#include <cstdlib>
#include <limits.h>
#include <getopt.h>
#include "RAN3.C"
#include "GASDEV_RAN3.C"

int main(int argc, char *argv[])
{
	float ndev=0., min=1., max=-1.;
	long int idum, nit=1L;
	int i=0, optcount=0;
	char c;
/*    	while((c=getopt(argc-1,argv+1,"n:"))!=-1){
		switch(c){
			case 'n':
				sscanf(optarg,"%ld",&nit);
				optcount+=1;
				break;
			default:
				break;
		}
	}
*/
	if((argc<=3+optcount))
	{
		fprintf(stderr, "\n\nusage: ndev <mean> <stdev> $RANDOM [niter] [min max]\n\n");
		exit(0);
	}
	sscanf(argv[3+optcount],"%ld", &idum);
	if(argc>4){
                sscanf(argv[4+optcount],"%ld", &nit);
	}
	if(argc>5){
                sscanf(argv[5+optcount],"%f", &min);
                sscanf(argv[6+optcount],"%f", &max);
		if(min>max){
			fprintf(stdout, "Minimum is greater than maximum!\n");
			fprintf(stderr, "\n\nusage: ndev <mean> <stdev> $RANDOM [niter] [min max]\n\n");
			exit(0);
		}
        }

	if(idum>0)
		idum*=-1;
	i=0;
	//fprintf(stdout, "%f %f %f\n", atof(argv[1+optcount]), atof(argv[2+optcount]),  atof(argv[3+optcount]));
	while (i<nit){
		ndev = (atof(argv[2+optcount])*gasdev(&idum))+atof(argv[1+optcount]);
		if(min<max){
			while ( ndev < min || ndev > max )
			{
				ndev = (atof(argv[2+optcount])*gasdev(&idum))+atof(argv[1+optcount]);
			}	
		}
		i++;
	}
	fprintf(stdout, "%f\n", ndev);
	return(0);
}
