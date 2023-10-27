#include "udf.h"
#include "unsteady.h"

DEFINE_PROFILE(inletBC_60C, t, nv)
{
    face_t f;

    
    if (CURRENT_TIME >= 0.0 && CURRENT_TIME < 6000.0)
    {
        begin_f_loop(f, t)
        {
            F_PROFILE(f, t, nv) = 0.08;
        }
        end_f_loop(f, t)
    }

    else if (CURRENT_TIME >= 6000.0)
    {
        begin_f_loop(f, t)
        {
            F_PROFILE(f, t, nv) = 0.0;
        }
    }
    
}