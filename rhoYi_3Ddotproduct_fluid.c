#include "udf.h"
#include "sg.h"
#include "mem.h"
#include "math.h"

DEFINE_PROFILE(int_bc_Yf_dotprod, t, nv)
{
    face_t f;
    real At[ND_ND], es[ND_ND], dr0[ND_ND], dG[ND_ND], dr1[ND_ND];
    real ds;
    real ds0, ds1;
    real l0, l1;
    real A_by_es;
    real A_by_es0, A_by_es1;
    real NV_VEC(es0), NV_VEC(es1);
    real S0, S1, h0, h1, gamma0, gamma1;
    real A[ND_ND];
    Thread* t0, * t1;
    real beta0, beta1;
    real alpha0, alpha1;
    real xf[ND_ND];
    real flux;
    int flag = 1;

    begin_f_loop(f, t)
    {
        cell_t c0 = F_C1(f, t);         /*  Solid side of the face  */
        cell_t c1 = F_C0(f, t);         /*  Fluid side of the face   */
        t0 = THREAD_T1(t);           /*  Solid side cell thread  */
        t1 = THREAD_T0(t);           /*  Fluid side cell thread    */
        gamma0 = C_DIFF_L(c0, t0, 0, 1);    /*  Solid side diffusion coefficient    */
        gamma1 = C_DIFF_EFF(c1, t1, 0);    /*  Fluid side diffusion coefficient  */
        real xc0[ND_ND];        /*  Solid cell centroid   */
        real xc1[ND_ND];        /*  Fluid cell centroid    */
        real phi0;
        real phi1;

        BOUNDARY_FACE_GEOMETRY(f, t, At, ds, es, A_by_es, dr0);      /* Check if A here is in the opposite of flux direction    */
        NV_D(A, =, (0 - At[0]), (0 - At[1]), (0 - At[2]));
        C_CENTROID(xc0, c0, t0);            /*  Solid cell centroid coordinates */
        C_CENTROID(xc1, c1, t1);            /*  Fluid cell centroid coordinates   */
        F_CENTROID(xf, f, t);

        ds0 = sqrt(pow((xf[0] - xc0[0]), 2) + pow((xf[1] - xc0[1]), 2) + pow((xf[2] - xc0[2]), 2));
        ds1 = sqrt(pow((xc1[0] - xf[0]), 2) + pow((xc1[1] - xf[1]), 2) + pow((xc1[2] - xf[2]), 2));
        NV_D(es0, =, (xf[0] - xc0[0]), (xf[1] - xc0[1]), (xf[2] - xc0[2]));     /*  Unit vector pointing from solid cell center to face */
        NV_D(es1, =, (-xc1[0] + xf[0]), (-xc1[1] + xf[1]), (-xc1[2] + xf[2]));  /*  Unit vector pointing from fluid cell center to face   */
        A_by_es0 = NV_DOT(A, A) / NV_DOT(A, es0);
        A_by_es1 = NV_DOT(A, A) / NV_DOT(A, es1);
        h0 = gamma0 * A_by_es0 / ds0;
        h1 = -gamma1 * A_by_es1 / ds1;

        phi0 = C_YI(c0, t0, 0);     /*  Cell value in solid side   */
        phi1 = C_YI(c1, t1, 0);     /*  Cell value in fluid side   */


        /*Message("Af = %.9f, %.9f, %.9f", A[0], A[1], A[2]);*/


        F_PROFILE(f, t, nv) =
            ((h0 * phi0 + h1 * phi1) / (h0 + h1));

        /*if (flag)
        {
            Message("Facevalue_fluid_side = %e\n;", F_PROFILE(f, t, nv));
        }
        flag = 0;*/

    }
    end_f_loop(f, t)
}