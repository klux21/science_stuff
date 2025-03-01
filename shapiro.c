#if 0
rm -f ./_shapiro*

if test "$(uname)" = "Linux"
then
   cc -Wall -O3 -o _shapiro ./shapiro.c -lm
else
   cc -Wall -O3 -o _shapiro ./shapiro.c
fi

./_shapiro
exit $?
#endif


/*****************************************************************************\
*                                                                             *
*  FILENAME:      shapiro.c                                                   *
*                                                                             *
* --------------------------------------------------------------------------- *
*                                                                             *
*  DESCRIPTION:   calculation of Shapiro delay according to the gravitational *
*                 potential of the Sun between Earth and Venus assuming the   *
*                 space is stretched according to the square of the Lorentz   *
*                 factor of the escape velocity                               *
*                                                                             *
* --------------------------------------------------------------------------- *
*                                                                             *
*  COPYRIGHT:     (c) 2024 Dipl.-Ing. Klaus Lux (Aachen, Germany)             *
*                                                                             *
* --------------------------------------------------------------------------- *
*                                                                             *
*  ORIGIN:        https://github/klux21/science_stuff                         *
*                                                                             *
* --------------------------------------------------------------------------- *
*                                                                             *
* Civil Usage Public License, Version 1.1, January 2024                       *
*                                                                             *
* Redistribution and use in source and binary forms, with or without          *
* modification, are permitted provided that the following conditions are met: *
*                                                                             *
* 1. Redistributions of source code must retain the above copyright           *
*    notice, this list of conditions, the explanation of terms                *
*    and the following disclaimer.                                            *
*                                                                             *
* 2. Redistributions in binary form must reproduce the above copyright        *
*    notice, this list of conditions and the following disclaimer in the      *
*    documentation or other materials provided with the distribution.         *
*                                                                             *
* 3. All modified files must carry prominent notices stating that the         *
*    files have been changed.                                                 *
*                                                                             *
* 4. The source code and binary forms and any derivative works are not        *
*    stored or executed in systems or devices which are designed or           *
*    intended to harm, to kill or to forcibly immobilize people.              *
*                                                                             *
* 5. The source code and binary forms and any derivative works are not        *
*    stored or executed in systems or devices which are intended to           *
*    monitor, to track, to change or to control the behavior, the             *
*    constitution, the location or the communication of any people or         *
*    their property without the explicit and prior agreement of those         *
*    people except those devices and systems are solely designed for          *
*    saving or protecting peoples life or health.                             *
*                                                                             *
* 6. The source code and binary forms and any derivative works are not        *
*    stored or executed in any systems or devices that are intended           *
*    for the production of any of the systems or devices that                 *
*    have been stated before except the ones for saving or protecting         *
*    peoples life or health only.                                             *
*                                                                             *
* The term 'systems' in all clauses shall include all types and combinations  *
* of physical, virtualized or simulated hardware and software and any kind    *
* of data storage.                                                            *
*                                                                             *
* The term 'devices' shall include any kind of local or non-local control     *
* system of the stated devices as part of that device als well. Any assembly  *
* of more than one device is one and the same device regarding this license.  *
*                                                                             *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   *
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     *
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  *
* POSSIBILITY OF SUCH DAMAGE.                                                 *
*                                                                             *
\*****************************************************************************/


#include <stdio.h>
#include <math.h>

#define CHECK_DIRECTION_DEPENDENCY 0  /* check direction dependency of the streching of the room
                                         that was propagated in German Wikipedia once (most likely wrong!) */

const double  year       = 31556952.0;     /* a year of 365.2425 days in seconds */
const double  day        = 86400.0;        /* a day in seconds */
const double  AE         = 149597870700.0; /* astromic unit in m */
const double  c          = 299792458.0;    /* speed of light m/s */
const double  G          = 6.6743e-11;     /* Newtonian const of gravitation m³/kg/s² */
const double  m_Sun      = 1.9884e30;      /* mass of sun in kg */
const double  r_Sun      = 696342000;      /* equatorial radius of sun */

const double  m_Earth    = 5.972168E24;    /* mass of Earth */
const double  r_Earth    = 6378137.0;      /* equatorial radius of Earth */
const double  r_Venus    = 6051800.0;      /* equatorial radius of Venus */
const double  rb_Mercury = 0.387098 * AE;  /* average radius of orbit of Mercury */
const double  rb_Pluto   = 39.482   * AE;  /* average radius of orbit of Pluto */
const double  rb_Venus   = 0.7233   * AE;  /* average radius of orbit of Venus */
const double  rb_Earth   = AE;             /* average radius of orbit of Earth */

const double  m_Andromeda  = 1.5e12  * m_Sun;
const double  rb_Andromeda = 2.5e6   * year * c;
//const double  m_Milkyway   = 1.15e12 * m_Sun;
const double  rb_Milkyway  = 2.67e4  * year * c;
const double  m_Milkyway   = 2.06e11 * m_Sun;
const double  rb_M87       = 53.5e6 * year * c;
const double  m_M87        = 6.0e12 * m_Sun;
const double  m_SGW        = 3.0*1.2e17 * m_Sun; /* estimated mass of Sloan Great Wall */
const double  rb_SGW       = 1.37e12 * year * c; /* estimated distance of Sloan Great Wall */
 
double sqr(double x)
{
    return(x*x);
}

int CalculateShapiroDelay()
{
   int iret = 0;

   double vb_Venus = sqrt(G * m_Sun / rb_Venus);            /* orbital velocity of Venus */
   double vb_Earth = sqrt(G * m_Sun / rb_Earth);            /* orbital velocity of Earth */
   double y_Venus  = 2.0 * M_PI * rb_Venus / vb_Venus;      /* calculated year of Venus */
   double y_Earth  = 2.0 * M_PI * rb_Earth / vb_Earth;      /* calculated year of Earth */
   double ys_Venus = 1.0 / (1.0 / y_Venus - 1.0 / y_Earth); /* synodic year of Venus (time it takes for Venus move around the sun from view of Earth) */
   double dist_Sun = r_Sun + r_Venus * 2.0;                 /* minimal distance from sun the radar may work; */
   double arc      = M_PI / 2.0 - asin(dist_Sun / rb_Venus) + M_PI / 2.0 - asin(dist_Sun / rb_Earth); /* angle between Earth and Venus from view of sun */
   double days     = 0.0;                                   /* days after visibility of Venus */

#if 1
   double v2_Earth = 0.0; /* -> ignore escape velocity from Earth */
#else
   /* enable this care about the difference in gravitational potential of the telescope on Earth and the potential on the way of light only */
   double v2_Earth = 2.0 * G * m_Sun / rb_Earth + 2.0 * G * m_Earth / r_Earth; /* square of escape velocity from Earth to care about */
#endif

   printf("Synodic year of Venus in calculation=%.3f years (%.4f days)\n", ys_Venus / year, ys_Venus / day);
   printf("Starting angle is %.4f° for passing the Sun at %.4f times the radius of the Sun\n", arc * 360.0 / 2.0 / M_PI, dist_Sun/r_Sun);
   
   /* calculation of Shapiro delay */
   while(arc >= 0.0)
   {
      double r        = 0.0; /* current distance from Earth or Venus */
      double sum      = 0.0; /* summary of additional distance */
#if CHECK_DIRECTION_DEPENDENCY
      double tan_sum  = 0.0; /* result if assuming lowered tangential speed that some people claim to exist */
      double diff_sum = 0.0; /* difference if assuming lower tangential speed that some people claim */
#endif
      double dist_max = 0.0; /* distance between one of the planets and the point of the connection line that the Sun is perpendicular to */

      double dist_Venus = sqrt(rb_Venus * rb_Venus + rb_Earth * rb_Earth - 2.0 * rb_Venus * rb_Earth * cos(arc)); /* distance between Earth and Venus to according the angle and the law of cosines */

      dist_Sun = sin(arc) / dist_Venus * rb_Earth * rb_Venus; /* minimal distance of the Sun to the ray according to the law of sines and that it just depends on the sinus
                                                                 of angle between the Sun and the other planet multiplied by the distance of the sun */

      printf("\n");
      printf("%.0f days after starting day when the angle between Earth and Venus is %.4f°\n", days, arc * 360.0 / 2.0 / M_PI);
      printf("distance of Venus is %.0f km\n", dist_Venus / 1000.0);
      printf("distance of passing the Sun is %.0f km (%.04f times the radius of Sun%s)\n",
              dist_Sun / 1000.0, dist_Sun / r_Sun,  arc >= (0.5 * M_PI) ? "" : ", behind Venus" );

      /* We know the distance of the sun that is perpendicular the direction to Venus and radius of the simplified circular orbits of the planets.
         Now we iterate the additional lengths from the point where that distance of the sun hits the path ouf our radar waves.
         Let's begin with the part the direction of the Venus. */

      if(arc <= (0.5 * M_PI))
         dist_max = 0.0; /* We need to skip this because Venus is closer to the Earth than our crossing point in this case. */
      else
         dist_max = sqrt(rb_Venus * rb_Venus - dist_Sun * dist_Sun);

      r = dist_max;
      while(r > r_Venus)
      {
         double cur_pos = dist_max - r;
         double dist    = sqrt (cur_pos * cur_pos + dist_Sun * dist_Sun); /* current distance of sun */
         double sqr_v2  = 2.0 * G * m_Sun / dist - v2_Earth; /* escape velocity at current position from sun */
         sum += 1000.0 * (sqr_v2 / (c*c - sqr_v2)); /* length increasement because of Lorentz factor */

#if CHECK_DIRECTION_DEPENDENCY
         /* dist_Sun / dist is the sine of the angle between the direction of movement according to the law of sines
            and the fact that dist is the site opposed the right angle which sine is 1. The sine of that angle multiplied
            by the length of that vector is the part of the vector that is perpendicular to the direction of Sun. */
         if(sqr_v2 > 0.0)
         {
            double tan_dist = 1000.0 * dist_Sun / dist * (sqrt(1.0 + sqr_v2 / (c*c - sqr_v2)) - 1.0);
            double vert_dist = 1000.0 * cur_pos / dist * (sqr_v2 / (c*c - sqr_v2));
            tan_sum += tan_dist;
            diff_sum += sqrt(sqr(tan_dist) + sqr(vert_dist));
         }
#endif

         if(r > 1000.0)
            r -= 1000.0;
         else
            r = 0.0;
      }

      /* Let's iterate the part towards Earth now. */
      dist_max = rb_Earth * rb_Earth - dist_Sun * dist_Sun;

      if(dist_max <= 0.0)
         dist_max = 0.0; /* ensure this */
      else
         dist_max = sqrt(dist_max);

      r = dist_max;

      if(arc < 0.5 * M_PI)
         r = dist_Venus - r_Venus; /* We have to iterate from current position of Venus only and to ignore the rest */

      while(r > r_Earth)
      {
         double cur_pos = dist_max - r;
         double dist    = sqrt (cur_pos * cur_pos + dist_Sun * dist_Sun); /* current distance of sun */
         double sqr_v2  = 2.0 * G * m_Sun / dist - v2_Earth; /* escape velocity at current position from sun */
         sum += 1000.0 * (sqr_v2 / (c*c - sqr_v2)); /* length increasement because of Lorentz factor */

#if CHECK_DIRECTION_DEPENDENCY
          /* dist_Sun / dist is the sine of the angle between the direction of movement according to the law of sines
             and the fact that dist is the site opposed the right angle which sine is 1. The sine of that angle multiplied
             by the length of that vector is the part of the vector that is perpendicular to the direction of Sun. */
         if(sqr_v2 > 0.0)
         {
            double tan_dist = 1000.0 * dist_Sun / dist * (sqrt(1.0 + sqr_v2 / (c*c - sqr_v2)) - 1.0);
            double vert_dist = 1000.0 * cur_pos / dist * (sqr_v2 / (c*c - sqr_v2));
            tan_sum += tan_dist;
            diff_sum += sqrt(sqr(tan_dist) + sqr(vert_dist));
         }
#endif
         if(r > 1000.0)
            r -= 1000.0;
         else
            r = 0.0;
      }

      sum *= 2.0;

      printf("summary additional distance=%.2f m resulting delay=%.2f us\n", sum, sum / c * 1000000.0);

#if CHECK_DIRECTION_DEPENDENCY
      diff_sum *= 2.0;
      printf("difference of distance=%.2f m resulting delay=%.2f us and summary delay=%.2f us\n",
             tan_sum, tan_sum / c * 1000000.0, diff_sum / c * 1000000.0);
#endif

      printf("\n");

      r = 2.0 * G * m_Andromeda / rb_Andromeda; /* square of escape velocity */
      sum = dist_Venus * r / (c*c - r); 
      printf("escape velocity from Andromeda galaxy is %.3f km/s causing an additional distance of %.0fm and a delay of %.6fus\n", sqrt(r) / 1000.0 , sum, sum / c * 1.0e6 );

      r = 2.0 * G * m_Milkyway / rb_Milkyway; /* square of escape velocity */
      sum = dist_Venus * r / (c*c - r); 
      printf("escape velocity from Milkyway galaxy is %.3f km/s causing an additional distance of %.0fm and a delay of %.6fus\n", sqrt(r) / 1000.0 , sum, sum / c * 1.0e6 );

      r = 2.0 * G * m_M87 / rb_M87; /* square of escape velocity */
      sum = dist_Venus * r / (c*c - r);
      printf("escape velocity from M87 galaxy is %.3f km/s causing an additional distance of %.0fm and a delay of %.6fus\n", sqrt(r) / 1000.0 , sum, sum / c * 1.0e6 );

      r = 2.0 * G * m_SGW / rb_SGW; /* square of escape velocity */
      sum = dist_Venus * r / (c*c - r);
      printf("escape velocity from Sloan Great Wall is %.3f km/s causing an additional distance of %.0fm and a delay of %.6fus\n", sqrt(r) / 1000.0 , sum, sum / c * 1.0e6 );

      if(arc == 0.0)
      { /* we are done and have finished our calculations now */
         arc = -1;
      }
      else
      {
         /* increase the number of days and the angle according to that time */
         days += 5.0;
         arc -= 2.0 * M_PI * 5.0 * day / ys_Venus;

         if(arc < 0.0)
         { /* Calculate the case that Venus is straight in the direction of the Sun finally */
            days += arc / (2.0 * M_PI) * ys_Venus / day;
            arc = 0.0;
         }
      }
   }

   iret = 1;
   return(iret);
}/* CalculateShapiroDelay() */


/* ------------------------------------------------------------------------- *\
   main function
\* ------------------------------------------------------------------------- */

int main(int argc, char * argv[])
{
   int iret = 0;

   if(!CalculateShapiroDelay())
      iret = 1;

   return (iret);
}/* main() */

/* ========================================================================= *\
   END OF FILE
\* ========================================================================= */
