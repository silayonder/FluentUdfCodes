#include "udf.h"
#include "sg.h"
#include "mem.h"
#include "math.h"


#define dpore 0.000000004   /*  Pore diameter for particle in m */
#define R 8.314 /*  Gas constant in J/molK    */
#define MW 0.046    /*  Molar mass of EtOH in kg/mol    */
#define MWA 46.0  /*    Molar mass of EtOH in g/mol */
#define MWB 28.0 /* Molar mass of N2 in g/mol    */
#define VA 50.36 /* Diffusion volume of EtOH */
#define VB 17.9 /*  Diffusion volume of N2    */
#define P 1.01325 /*    1 atm   */
#define pi 3.14159


DEFINE_DIFFUSIVITY(diff_coeff, c, thread, i)
{
    real diff;
    real sphereID = 13;   /*Solid ID for particle cell zone*/
    real fluidID = 12;    /*Fluid ID for fluid cell zone*/
    
    if (THREAD_ID(thread) == sphereID)
    {
        diff = (dpore/3.0) * sqrt((8.0*R*C_T(c, thread))/(pi*MW));
    }   
        
    else 
    {
        diff = ((0.001*pow(C_T(c, thread), 1.75) * sqrt(1.0/MWA + 1.0/MWB)) / (P*(pow(VA, 1.0/3.0) + pow(VB, 1.0/3.0))*(pow(VA, 1.0/3.0) + pow(VB, 1.0/3.0)))) * 0.0001;
    }
        
        return diff;
}
