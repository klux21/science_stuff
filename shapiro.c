#if 0
rm -f ./_shapiro*
cc -Wall -O3 -o _shapiro ./shapiro.c
./_shapiro
exit $?
#endif


/*****************************************************************************\
*                                                                             *
*  FILENAME:      shapiro.c                                                   *
*                                                                             *
* --------------------------------------------------------------------------- *
*                                                                             *
*  DESCRIPTION:   calculation of Shapiro delay according to the difference of *
*                 gravitational potential of a telescope on Earth and Venus   *
*                 assuming the space is stretched according to the square of  *
*                 the Lorentz factor of the speed of a free fall between      *
*                 between the potentials                                      *
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

const double AE         = 149597870700.0;/* astromic unit in m */
const double c          = 299792458.0;   /* speed of light m/s */
const double G          = 6.6743e-11;    /* Newtonian const of gravitation m³/kg/s² */
const double m_Sun      = 1.9884e30;     /* mass of sun in kg */
const double r_Sun      = 696342000;     /* equatorial radius of sun */

const double m_Earth    = 5.972168E24;   /* mass of Earth */
const double r_Earth    = 6378137.0;     /* equatorial radius of Earth */
const double r_Venus    = 6051800.0;     /* equatorial radius of Venus */
const double rb_Mercury = 0.387098 * AE; /* average radius of orbit of Mercury */
const double rb_Pluto   = 39.482   * AE; /* average radius of orbit of Pluto */
const double rb_Venus   = 0.7233   * AE; /* average radius of orbit of Venus */
const double rb_Earth   = AE;            /* average radius of orbit of Earth */


/* ------------------------------------------------------------------------- *\
   main function
\* ------------------------------------------------------------------------- */

int main(int argc, char * argv[])
{
   int iret = 0;

   double dist_Sun = r_Sun + r_Venus;
   double v2_Earth = 2.0 * G * m_Sun / rb_Earth + 2.0 * G * m_Earth / r_Earth;
   double r = 0.0;
   
   /* calculation of Shapiro delay */
   while(dist_Sun <= rb_Venus)
   {
       double dist_max = rb_Venus * rb_Venus - dist_Sun * dist_Sun;
       double sum = 0.0;
       double diff_sum = 0.0; /* difference if using lowered tanential speed */

       if(dist_max <= 0.0)
          dist_max = 0.0; /* ensure this */
       else
          dist_max = sqrt(dist_max);

       r = dist_max;
       while(r > 0.0)
       {
          double cur_pos = dist_max - r;
          double dist    = sqrt (cur_pos * cur_pos + dist_Sun * dist_Sun); /* current distance of sun */
          double sqr_v2  = 2.0 * G * m_Sun / dist - v2_Earth; /* escape velocity at current position from sun */
          sum += 1000.0 * (sqr_v2 / (c*c - sqr_v2)); /* length increasement because of Lorentz factor */

          if(sqr_v2 > 0)
             diff_sum += 1000.0 * dist_Sun / dist * (sqrt(1.0 + sqr_v2 / (c*c - sqr_v2)) - 1.0);

          if(r > 1000.0)
             r -= 1000.0;
          else
             r = 0.0;
       }

       dist_max = rb_Earth * rb_Earth - dist_Sun * dist_Sun;

       if(dist_max <= 0.0)
          dist_max = 0.0; /* ensure this */
       else
          dist_max = sqrt(dist_max);

       r = dist_max;
       while(r > 0.0)
       {
          double cur_pos = dist_max - r;
          double dist    = sqrt (cur_pos * cur_pos + dist_Sun * dist_Sun); /* current distance of sun */
          double sqr_v2  = 2.0 * G * m_Sun / dist - v2_Earth; /* escape velocity at current position from sun */
          sum += 1000.0 * (sqr_v2 / (c*c - sqr_v2)); /* length increasement because of Lorentz factor */

          if(sqr_v2 > 0)
             diff_sum += 1000.0 * dist_Sun / dist * (sqrt(1.0 + sqr_v2 / (c*c - sqr_v2)) - 1.0);

          if(r > 1000.0)
             r -= 1000.0;
          else
             r = 0.0;
       }

       sum *= 2.0;
       diff_sum *= 2.0;

       printf("\n");
       printf("summary additional distance=%.2f m resulting delay=%.8f s\n", sum, sum / c);
       printf("difference of distance=%.2f m resulting delay=%.8f s resulting summary delay=%.8f s\n", diff_sum, diff_sum / c, (sum - diff_sum) / c );

       if (dist_Sun == rb_Venus)
         dist_Sun *= 2.0;
       else
         dist_Sun = rb_Venus;
   }

   return (iret);
}/* main() */

/* ========================================================================= *\
   END OF FILE
\* ========================================================================= */
