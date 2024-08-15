#include "udf.h"
#include "sg.h"
#include "mem.h"
#include "math.h"
#include "stdio.h"

#define sphereID 37796
#define W0 0.3955
#define D 0.0006049
#define n 1.156
#define A 5.24677
#define B 1598.673
#define C -46.424
#define R 0.00008314 /*Gas constant in m3bar/molK*/
#define MW 0.046    /*Molar mass of EtOH in kg/mol*/

DEFINE_EXECUTE_AT_END(execute_sourceIntegral)
{
	real sum_source = 0.0;
	real ksav;
	real W;
	real Ps;
	real Pi;
	real w;
	real source;

#if !RP_HOST /*	Get the variables that will be used in computation and are not stored on the host */
	Thread* t;
	cell_t c;
	Domain* d;
	d = Get_Domain(1);
#endif

#if !RP_HOST

	thread_loop_c(t, d)
	{
		if (THREAD_ID(t) == sphereID)
		{
			begin_c_loop_int(c, t)		/*	Integrate over interior cells of the partitioned mesh (not to compute some nodes double)	*/
			{
				Ps = pow(10.0, (A - (B / (C_T(c, t) + C))));    /*Calculate vapor pressure by Antoine in bar*/
				Pi = C_R(c, t) * C_YI(c, t, 0) * C_T(c, t) * R / MW;
				W = W0 * exp(-D * pow((C_T(c, t) * log(Ps / (Pi + SMALL))), n));  /*Adsorption capacity in g/g*/
				ksav = 0.58 * exp(-1970.0 / C_T(c, t));
				w = C_UDSI(c, t, 0);    /*Amount adsorbed in g/g*/

				if (Pi == 0.0)
				{
					source = 0.0;
				}

				else
				{
					source = ksav * (W - w);
				}

				sum_source += source * C_VOLUME(c, t);

			}
			end_c_loop_int(c, t)

		}
	}

# if RP_NODE	/*	Compute the sum of the integral across all nodes	*/

	sum_source = PRF_GRSUM1(sum_source);

# endif

#endif

	node_to_host_real_1(sum_source);	/*	Send the total integral to host for printing	*/

#if !RP_NODE
	real time = RP_Get_Real("flow-time");
	FILE* fp = NULL;
	char filename[] = "DenemeSourceParallel_Packed_YeniDeneme.txt";
	fp = fopen(filename, "a");
	fprintf(fp, "%g, %.16f\n", time, sum_source);
	fflush(fp);
	fclose(fp);
#endif

}