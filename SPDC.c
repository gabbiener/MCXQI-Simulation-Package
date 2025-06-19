/*
* Spdc
NOTE: NO phase change upon reflection included yet!!!!
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "ray.h"
#include "rng-double.h"
#define PI 3.14159265358979323846264338327


enum {
        GAUSSIAN,
        UNIFORM,
        NONE
};

enum {
	R_SEED,
	R_DOUBLE,
	R_FREE
};

double as2a(double AlphaS, double a)
{
  return 2*acos(a*sin(AlphaS)/sqrt(pow(a,2)-2*a*cos(AlphaS)+1));
}

double as2ai(double AlphaS, double a)
{
  //fprintf(stderr, "AlphaI: %le, AlphaS: %le\n", acos((2*a-(1+pow(a,2))*cos(AlphaS))/(1+pow(a,2)-2*a*cos(AlphaS))), AlphaS);
  return acos((2*a-(1+pow(a,2))*cos(AlphaS))/(1+pow(a,2)-2*a*cos(AlphaS)));
}

double Calc_AlphaS_Min(double maxAlphaS, double a)
{
  double Alpha = as2a(maxAlphaS, a);
  return Alpha - maxAlphaS;
}

double I_preCDF(double AlphaS, double a)
{
  return 1/sqrt(1-pow(a,2))*atan(sqrt((1+a)/(1-a))*tan(AlphaS/2))-a*sin(AlphaS)/2/(1-a*cos(AlphaS));
}

double I_DRV(double AlphaS, double a)
{
  return 1/(1-a*cos(AlphaS))-(1-pow(a,2))/2/pow((1-a*cos(AlphaS)),2);
}

double Newton_Raphson(double u, double a, double minAlphaS, double maxAlphaS, double Eps, int nIter)
{
  double AlphaS      = acos(a);
  double AlphaS_Prev = AlphaS;
  double I_minAlphaS = I_preCDF(minAlphaS, a);
  double I_maxAlphaS = I_preCDF(maxAlphaS, a);
  double Norm_Factor = I_maxAlphaS - I_minAlphaS;
  double I_AlphaS    = 0;
  //double F           = 0;
  //double f           = 0;
  //fprintf(stderr, "minAlphaS: %le, minAlphaS: %le \n", minAlphaS, maxAlphaS);
  for (int j = 1; j <= nIter; j++)
    {
      if (j==nIter)
	fprintf(stderr,"Warning: loop counter has exceeded the maximum limit!\n");
      
      I_AlphaS = I_preCDF(AlphaS, a);
      AlphaS_Prev = AlphaS;
      AlphaS = AlphaS - (I_AlphaS-I_minAlphaS-u*Norm_Factor)/I_DRV(AlphaS, a);
      if (fabs(AlphaS_Prev - AlphaS) < Eps)
	{
	  break;
	}
    }
    
  return AlphaS;
    
}

double* rand_double(int c, int n)
{
/* 
   An encapsulation object that returns a pointer to an array of double precision
   random numbers.  The idea is to decouple the random number generator from the
   rest of the code, making a change of random number generators easier.
   if c = R_SEED, then n = seed value
   if c = R_DOUBLE, then n = # of returned doubles
   if c = R_FREE, then n = void
*/	
	static int length = 1024; //default length
	static int index = -1;
	static double *rand_array = NULL;
	
	switch (c) {

		case R_SEED:
		ranf_start(n);
		index = 0;
		return NULL;
		break;

		case R_DOUBLE:
		if (index < 0) //initialize 
			rand_double(R_SEED, time(NULL));
		if (rand_array == NULL) 
			rand_array = malloc(length * sizeof(double));
			// NOTE: should really check return value from malloc
		if (n >= length) {
			free(rand_array);
			length = n;
			rand_array = malloc(length * sizeof(double));
			index = 0;
		}			
		if ((index == 0) || (n + index >= length)) {
			ranf_array(rand_array,length);
			index = 0;
		}
		index = index + n;		
		return rand_array + index - n; 
		//ok I could add a variable to make math easier
		break;

		case R_FREE:
		free(rand_array);
		return NULL;
		
		default:
		return NULL;
	}
	return NULL;	
}


static double gauss_rand()
{

        //creats 0 mean & unit sigma gaussian random numbers    
        //using Box-Muller transfomre (there are better methods)
        double dist, mag;
        double *var;
        do {
                var = rand_double(R_DOUBLE,2);
                var[0] = 2.0*var[0]-1.0;
                var[1] = 2.0*var[1]-1.0;
                dist = (var[0]*var[0]) + (var[1]*var[1]);
        } while (dist >= 1.0 || dist == 0.0);
        mag = sqrt(-2.0*log(dist)/dist);
        dist = mag*var[0];
        return dist;
}

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);

printf("Creates a parallelogram crystal and propagates all rays to the crystal.\n"
"Some rays in crystal aperture decay into two smaller energy rays, and it sets all rays outside\n"
"crystal aperture to 0 intensity. The crystal is defined by the parallelogram\n"
"with vertices p0, p1, p2, and (p0 + p1 + p2).\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help                     Display this help message.\n"
"  -i, --input=<file>             Input filename. Default: stdin.\n"
"  -o, --output=<file>            Output filename. Default: stdout.\n"
"      --P0='[x,y,z]'             First vertex (default = [0,0,0].\n"
"      --P1='[x,y,z]'             Second vertex (default = [1,0,0].\n"
"      --p2='[x,y,z]'             Third vertex (default = [0,1,0].\n"
"  -m  --miller='[mh,mk,ml]'      Miller index (default = [1,1,1]. \n"
"  -l  --lattice=<value>          Lattice Constant (default = 3.57e-10). \n"
"  -p  --arange='[amin,amax]'     specified Alpha range (default = [0,PI].\n"
"  -r                             RockinCurve.\n"
"  --seed=<num>                   Seed for random number generator (default is current time)\n"
"  -z, --zero                     Enables propagation of 0 intensity rays\n"
"  -f, --RandomRotate             Switch K1 and k2 randomly\n"
"                                 By default not propaged.\n" 
"\n");
}

int main(int argc, char *argv[])
{

	/* initial variable */
	//file stuff
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	// input options stuff
	int c;
	double Eps = 1e-9;
	int nIter = 1000;
	int seed = time(NULL);
	int useZero = 0;
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);
	
	// Miller indices and lattice constant for the crystal
	int mh = 1; // Note: these are fixed, change later
	int mk = 1;
	int ml = 1;
	double latticeConstant = 3.57e-10;
	double minAlphaS = 0;
	double maxAlphaS = PI;
	double sigmarock = 0;
	double angRange = 2*PI;
	double angBias = 0;
	int RandomRotate = 0;
	// math stuff
	Vector B1, B2, B1unitV, B2unitV, Normal; //basis vectors & normal 
	double B1mag, B2mag;
	FILE *Alpha_file;


	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"P0",                1, NULL,               'a'},
		{"P1",                1, NULL,               'b'},
		{"P2",                1, NULL,               'c'},
		{"miller",            1, NULL,               'm'},
		{"lattice",           1, NULL,               'l'},
		{"arange",            1, NULL,               'p'},
		{"sigmarock",         1, NULL,               'r'},
		{"angRange",          1, NULL,               'g'},
		{"angBias",           1, NULL,               'd'},
		{"seed",              1, NULL,               's'},
		{"RandomRotate",      0, NULL,		     'f'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:m:l:p:r:g:d:s:z",
                longopts, NULL)) != -1)
        {

                switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			inFileName = strdup(optarg);
			break;

			case 'o' :
			outFileName = strdup(optarg);
			break;

			case 'a' :
			sscanf(optarg,"[%le,%le,%le]", &P0.x, &P0.y, &P0.z);
			break;

			case 'b' :
			sscanf(optarg,"[%le,%le,%le]", &P1.x, &P1.y, &P1.z);
			break;

			case 'c' :
			sscanf(optarg,"[%le,%le,%le]", &P2.x, &P2.y, &P2.z);
			break;
			
			case 'm':
   			sscanf(optarg, "[%d,%d,%d]", &mh, &mk, &ml);
    			break;

			case 'l':
                        latticeConstant = strtod(optarg,NULL);
                        break;
			
			case 'p' :
                        sscanf(optarg,"%le", &maxAlphaS);
                        break;

			case 'r':
                        sigmarock = strtod(optarg,NULL);
                        break;	

			case 'g':
                        angRange = strtod(optarg,NULL)*PI/180;
                        break;

			case 'd':
                        angBias = strtod(optarg,NULL)*PI/180;
                        break;
			
			case 'f':
                        RandomRotate = 1;
                        break;

			case 's' :
			seed = atoi(optarg);
			break;

			case 'z' :
			useZero = 1;
			break;

			default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;

		}
	}

	/* set random seed */
	rand_double(R_SEED, seed);
	
	//calculate cross product sanity check
	Normal = cross(sub(P1,P0),sub(P2,P0));
	
	if (magnitude(Normal) == 0.0) {  //see if points are colinear
		fprintf(stderr,"Error points are colinear!\n");
		return 1;
	}
	Normal = unit(Normal);

	//fprintf(stderr,"Normal: %g %g %g\n", Normal.x, Normal.y, Normal.z); //TEST

	/* open files if necessary */
	if ( inFileName != NULL ) {
		in = fopen(inFileName, "rb");
		if ( in == NULL ) {
			fprintf(stderr,"Failed to open input file: %s!\n", inFileName);
			return 1;
		}
		free(inFileName);
	}

	if ( outFileName != NULL ) {
		out = fopen(outFileName, "wb");
		if ( out == NULL ) {
			fprintf(stderr,"Failed to open output file: %s!\n", outFileName);
			return 1;
		}
		free(outFileName);
	}

	/* calculate unit basis vectors */
	B1 = sub(P1,P0);
	B1unitV = unit(B1);
	B1mag = magnitude(B1);
	B2 = sub(P2,P0);
	B2unitV = unit(B2);
	B2mag = magnitude(B2);

	/* open a file for recording p_values */	
	Alpha_file = fopen("Alpha_Values.txt", "w");
	if (Alpha_file == NULL) {
        printf("File does not exist!");
        return 1;
   	}		
    
    
	/* loop over vectors */
	while (!feof(in)) {				
		double distance, m, n, point_B1, point_B2, B1_B2;		
		Ray myRay;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		Vector point;
		
		myRay.v = unit(myRay.v); //make the direction a unit vector
		
		/* check if using 0 or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* check ray is not parallel to crystal */
		if (dot(myRay.v,Normal) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " crystal surface!\n Skipping ray, but copying it to output.\n", myRay.ray);
			write_ray(out, myRay);
			continue;
		} 

		//caluclate distance from the location point in myRay to 
		//aperture plane along pR_DOUBLEropagation direction
		distance = dot(sub(P0,myRay.pos),Normal)/dot(myRay.v,Normal);
		if (distance < 0.0) {
			fprintf(stderr,"Warning Ray %lu is furthern then the crystal's surface"
			 " back propagating to crystal.\n", myRay.ray);
		}
		myRay.pos = add(myRay.pos,mult(myRay.v,distance));
		myRay.l = myRay.l + distance;

		//calculate where myRay.pos is with respect to P0 & basis vectors
		//B1=(P1-P0), B2=(P2-P0).
		point = sub(myRay.pos, P0);
		point_B1 = dot(point,B1unitV);
		point_B2 = dot(point,B2unitV);
		B1_B2    = dot(B1unitV,B2unitV);		
		m = (point_B1-(point_B2*B1_B2))/(1-(B1_B2*B1_B2)); //Basis 1 length
		n = point_B2-(m*B1_B2); //Basis 2 length
		m = m/B1mag; // was using unit lenght basis vectors to make math eaier
		n = n/B2mag; // convert back to regular size vector
		
		if ( m >= 0.0 && m <= 1.0 && n >= 0.0 && n <= 1.0 ) {
		  //in crystal
		  double MagG, Magkp, MagkG, *d_rand, lambda1, lambda2, ALPHA;
		  Vector kp, G, u, Gprime, kG, kGunitV, k1prime, k2prime, k1, k2;
		  Ray myRay1, myRay2;
			Magkp = 2*PI/myRay.w;
			kp = mult(myRay.v, Magkp);
			MagG = (2*PI*sqrt(pow(mh,2)+pow(mk,2)+pow(ml,2))) / latticeConstant;
			G = mult(Normal,MagG);
      			G = unit(G);
			u = unit(sub(P1,P0));
			d_rand =  rand_double(R_DOUBLE,4);
			double drock_radian = sigmarock*PI/180;
			double cos_alpha = cos(PI*d_rand[0]);
			double sin_alpha = sin(PI*d_rand[0]);
			double gr = gauss_rand();
			double cos_beta = cos(drock_radian*gr);
			double sin_beta = sin(drock_radian*gr);
			Gprime.x = cos_beta*G.x + sin_beta*cos_alpha*(u.x*G.z-u.z*G.y) + sin_beta*sin_alpha*u.x;
			Gprime.y = cos_beta*G.y + sin_beta*cos_alpha*(u.z*G.x-u.x*G.z) + sin_beta*sin_alpha*u.y;
			Gprime.z = cos_beta*G.z + sin_beta*cos_alpha*(u.x*G.y-u.y*G.x) + sin_beta*sin_alpha*u.z;
			Gprime = mult(Gprime,MagG);
			kG = add(kp,Gprime);
			MagkG = magnitude(kG);
			double a = MagkG/Magkp;
			//fprintf(stderr, "Magkp, MagkG, a: %le, %le, %le\n", Magkp, MagkG, a);
		
			if ( MagkG > Magkp ) {			
				myRay.i = 1.0; 
				write_ray(out, myRay);
				// Remark outside SPDC condition

			} else {
			
			
				//Rotate kG parallel to z axis using rotation matrix R
				//NkG.x = (pow(kG.y,2)+pow(kG.x,2)*kG.z) * kG.x /(pow(kG.x,2)+pow(kG.y,2)) + kG.x*kG.y*(-1+kG.z)/(pow(kG.x,2)+pow(kG.y,2))*kG.y - kG.x * kG.z;
				//NkG.y = (kG.x*kG.y*(-1+kG.z))/(pow(kG.x,2)+pow(kG.y,2))*kG.x +  (pow(kG.x,2)+pow(kG.y,2)*kG.z) * kG.y/(pow(kG.x,2)+pow(kG.y,2)) - kG.y * kG.z;
				//NkG.z = kG.x*kG.x + kG.y*kG.y + kG.z*kG.z;
				//NkG.x = 0;
				//NkG.y = 0;
				//NkG.z = magnitude(kG);
								

				// Here we need to rotate k1 and k2 in the prime system about the z-axis randomly by 2Pi
				double cos_theta = cos(2*PI*d_rand[1]); 
				double sin_theta = sin(2*PI*d_rand[1]);
				
				//double theta = angRange*(d_rand[1]-0.5) + angBias; //angRange and angBias give us freedom to filtere illumination to a specific angle of cone if it desired; we flip signal and idler randomly by RandomRotate=1 
				//if (d_rand[2] > 0.5) {
				//	theta += RandomRotate*PI;
				//}
				//double cos_theta = cos(theta);   
				//double sin_theta = sin(theta);   
				double xr  = d_rand[2];
				//int Flip  = d_rand[3] > 0.5;

				
				maxAlphaS = fmin(maxAlphaS, PI);
				minAlphaS = Calc_AlphaS_Min(maxAlphaS, a);
				double AlphaS = Newton_Raphson(xr, a, minAlphaS, maxAlphaS, Eps, nIter);
				double AlphaI = as2ai(AlphaS, a);

			        //if (d_rand[1] < 0.5)
				//{
				//  fprintf(Alpha_file,"%le, %le, %le\n",AlphaS, AlphaI, AlphaS+AlphaI);
				//} else
				//{
				//  fprintf(Alpha_file,"%le, %le, %le\n",AlphaI, AlphaS, AlphaS+AlphaI);
				//}

				fprintf(Alpha_file,"%le, %le, %le\n",AlphaS, AlphaI, AlphaS+AlphaI);
				
				double h = MagkG*sin(AlphaI)*sin(AlphaS)/sin(AlphaS+AlphaI);

				k1prime.x = - sin_theta * h;
   				k1prime.y =   cos_theta * h;
	    			k1prime.z = MagkG*sin(AlphaI)*cos(AlphaS)/sin(AlphaS+AlphaI);
					
				k2prime.x = - sin_theta * (-h);
				k2prime.y =   cos_theta * (-h);
				k2prime.z = MagkG*sin(AlphaS)*cos(AlphaI)/sin(AlphaS+AlphaI);

				//return rotated k1 and k2 vectors using transpose of R matrix
				
				kGunitV = unit(kG);					

				k1.x = (pow(kGunitV.y,2)+pow(kGunitV.x,2)*kGunitV.z) /(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k1prime.x + (kGunitV.x*kGunitV.y*(-1+kGunitV.z))/(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k1prime.y + kGunitV.x*k1prime.z;
				k1.y = (kGunitV.x*kGunitV.y*(-1+kGunitV.z))/(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k1prime.x + (pow(kGunitV.x,2)+pow(kGunitV.y,2)*kGunitV.z) /(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k1prime.y + kGunitV.y*k1prime.z;
				k1.z = - kGunitV.x*k1prime.x - kGunitV.y*k1prime.y + kGunitV.z*k1prime.z;
				
				k2.x = (pow(kGunitV.y,2)+pow(kGunitV.x,2)*kGunitV.z) /(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k2prime.x + (kGunitV.x*kGunitV.y*(-1+kGunitV.z))/(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k2prime.y + kGunitV.x*k2prime.z;
                                k2.y = (kGunitV.x*kGunitV.y*(-1+kGunitV.z))/(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k2prime.x + (pow(kGunitV.x,2)+pow(kGunitV.y,2)*kGunitV.z) /(pow(kGunitV.x,2)+pow(kGunitV.y,2))*k2prime.y + kGunitV.y*k2prime.z;
                                k2.z = - kGunitV.x*k2prime.x - kGunitV.y*k2prime.y + kGunitV.z*k2prime.z;

					
				lambda1 = 2*PI/magnitude(k1);
				lambda2 = 2*PI/magnitude(k2);
				
				k1 = unit(k1);
				k2 = unit(k2);

				myRay1 = myRay;
 				myRay1.v.x = k1.x;
				myRay1.v.y = k1.y;
				myRay1.v.z = k1.z;
				myRay1.w = lambda1;

				myRay2 = myRay;
				myRay2.v.x = k2.x;
				myRay2.v.y = k2.y;
				myRay2.v.z = k2.z;
				myRay2.ray = myRay.ray + 1;
				myRay2.w = lambda2;

				write_ray(out, myRay1);
				write_ray(out, myRay2);
				

			}
	
                } else {
			myRay.i = 0.0; // ray is outside crystal aperture
			write_ray(out, myRay);
		}

		//print ray
		//PS write_ray(out, myRay);
	}

	fclose(in);
	fclose(out);
	fclose(Alpha_file);
	
	rand_double(R_FREE,0);

	return 0;
}
