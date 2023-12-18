#include "udf.h"
#include "sg.h"
#include "mem.h"
#include "math.h"

#define W0 0.3955
#define D 0.0006049
#define n 1.156
#define A 5.24677
#define B 1598.673
#define C -46.424
#define R 0.00008314 /*Gas constant in m3bar/molK*/
#define MW 0.046    /*Molar mass of EtOH in kg/mol*/
#define RHOS 750.0 /*Solid density of particle in kg/m3*/
#define EPS 0.35    /*Particle porosity*/

DEFINE_SOURCE(c_source, c, t, dS, eqn)
{
    real ksav;
    real W;
    real Ps;
    real Pi;
    real w;
    real source;
    real dwdt;

    Ps = pow(10.0, (A - (B / (C_T(c, t) + C))));    /*Calculate vapor pressure by Antoine in bar*/
    Pi = C_R(c, t) * R * C_YI(c, t, 0) * C_T(c, t) / MW;
    W = W0 * exp(-D * pow((C_T(c, t) * log(Ps / (Pi + SMALL))), n));  /*Adsorption capacity in g/g*/
    ksav = 0.58 * exp(-1970.0 / C_T(c, t));
    w = C_UDSI(c, t, 0);    /*Amount adsorbed in g/g*/
    dwdt = ksav * (W - w);

    if (Pi == 0)
    {
        source = 0;
    }

    else {
        source = -(1 - EPS) * RHOS * dwdt / (EPS);
    }


    return source;
}