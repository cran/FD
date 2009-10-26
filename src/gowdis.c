/* Etienne Laliberte, October 26, 2009 
 * etiennelaliberte[AT]gmail.com */

#include <R.h>
#include <math.h>

static int dequal(double d1, double d2)
{
  double diff = fabs(d1 - d2);
  return (diff < 0.00000001);
}

void gowdis(double *x, double *win, int *type, int *nin, int *pin, double *range2, int *podin, double *tmax, double *tmin, double *res)
{
   int i, j, k, n = *nin, p = *pin, pod = *podin, l;
   double w, wsum, psim, tmp, ti, tj;
   
   for (i = 0; i < n-1; i++)
   {
     for (j = i+1; j < n; j++)
     {
          wsum = 0.0; tmp = 0.0;
          for (k = 0; k < p; k++)
          {
             psim = 0.0; w = 0.0;
             if (R_FINITE(x[i + n * k]) && R_FINITE(x[j + n * k]) )
             {
               w = win[k];
               switch (type[k])
               {
                 case 1:
                 {
                   psim = w * (1 - (fabs(x[i + n * k] - x[j + n * k])) / range2[k]);
                   break;
                 }
                 case 2:
                 {
                   switch (pod)
                   {
                     case 1:
                     {
                       if (dequal(x[i + n * k], x[j + n * k])) 
                       {
                         psim = w * 1.0;
                       }
                       else
                       {
                          ti = 0.0; tj = 0.0;
                          for (l = 0; l < n; l++){
                              if (R_FINITE(x[l + n * k]) && dequal(x[l + n * k], x[i + n * k])) ti = ti + 1.0;
                              if (R_FINITE(x[l + n * k]) && dequal(x[l + n * k], x[j + n * k])) tj = tj + 1.0;
                           }
                           psim = w * (1- ((fabs(x[i + n * k] - x[j + n * k]) - (ti-1)/2 - (tj-1)/2  ) / (range2[k] - (tmax[k]-1)/2 - (tmin[k]-1)/2 ) ) ) ;
                       }
                       break;
                     }
                    case 2:
                    {
                      psim = w * (1 - (fabs(x[i + n * k] - x[j + n * k])) / range2[k]);
                      break;
                    }
                  }
                  break;
                 }
                 case 3:
                 case 4:
                 {
                  if (dequal(x[i + n * k], x[j + n * k])) psim = w * 1.0;
                   break;
                 }
                 case 5:
                 {
                  if (dequal(x[i + n * k], 0) && dequal(x[j + n * k], 0))
                  {
                     w = 0.0;
                     psim = 0.0;
                  }
                  else
                  {
                    if (dequal(x[i + n * k], x[j + n * k])) psim = w * 1.0;
                  }
                   break;
                 }
                 default:
                 {
                   psim = 0.0;
                   break;
                 }
               }
            }
            tmp += psim;
            wsum += w;
        }
     if (wsum == 0)  *res++ = NA_REAL;
     else *res++ = 1 - (tmp / wsum);
   }
  }
}

