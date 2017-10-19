/*
The MIT License (MIT)

Copyright (c) 2003-2017 Syoyo Fujita and many contributors.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#ifndef TINYSHROT_H_
#define TINYSHROT_H_

#include <malloc.h>
#include <math.h>

/*
 * Builds rotation matrix for spherical harmoinics.
 * `a`, `b` and `c` parameter are based on 
 * "Evaluation of the rotation matrices in the basis of real spherical harmonics", 1997
 *
 * TODO(syoyo): Require more tests and mathematical verification.
 * NOTE(syoyo): Not numerically stable for larger l(l <= 4 should work well)
 *
 * `ret` must be memory allocated like this:
 *
 * double **mat;
 * 
 * mat = (double **)malloc(sizeof(double *) * (l + 1));
 * 
 * for (i = 0; i <= l; i++) {
 * 
 *     mat[i] = (double *)malloc(sizeof(double) * (2 * l + 1) * (2 * l + 1));
 * 
 * }
 *
 * tinysh_rotation(mat, l, a, b, r);
 * 
 */

void tinysh_rotation(double **ret, unsigned int l, double a, double b, double r);

#endif /* TINYSHROT_H_ */

#ifdef TINYSHROT_IMPLEMENTATION
#ifndef TINYSHROT_IMPLEMENTATION_DEIFNED
#define TINYSHROT_IMPLEMENTATION_DEIFNED
#
static void buildcoeff(double *ret, double *c, int l, double a, double r);
static double phi(int m, double p);
static double sign(int m);
static int map(int l, int m, int md);
static void fill_d(double *m, int l);

int map(int l, int m, int md) {
  /* map subscript (m, m') to array index.
   * mapping is
   * f(m, m') = (m + l) * (2 * l + 1) + (m' + l)
   */
  return (m + l) * (2 * l + 1) + (md + l);
}

double sign(int m) {
  if (m >= 0)
    return 1.0;

  return -1.0;
}

double phi(int m, double p) {
  double val = 0.0;

  if (m > 0) {
    val = sqrt(2.0) * cos((double)(m * p));
  } else if (m < 0) {
    val = sqrt(2.0) * sin((double)(-m * p));
  } else {
    val = 1.0;
  }

  return val;
}

void tinysh_rotation(double **ret, unsigned int l, double a, double b, double r) {
  int i;
  int id;
  int mid;
  int mid1, mid2;
  int m, md;
  double **d;
  double c1, c2, c3;

  d = (double **)malloc(sizeof(double *) * (l + 1));

  for (i = 0; i <= l; i++) {
    d[i] = (double *)malloc(sizeof(double) * (2 * l + 1) * (2 * l + 1));
  }

  d[0][0] = 1.0;
  ret[0][0] = 1.0;

  id = map(1, 0, 0);
  d[1][id] = cos(b); /* d1_(0, 0)	*/

  id = map(1, 1, -1);
  d[1][id] = sin(b * 0.5) * sin(b * 0.5); /* d1_(1,-1)	*/

  id = map(1, 1, 0);
  d[1][id] = (-1.0 / sqrt(2.0)) * sin(b); /* d1_(1, 0)	*/

  id = map(1, 1, 1);
  d[1][id] = cos(b * 0.5) * cos(b * 0.5); /* d1_(1, 1)	*/

  /* d1_(0,-1) = d1_(1, 0) */
  id = map(1, 0, -1);
  d[1][id] = d[1][map(1, 1, 0)];

  /* d1_(0,1) = -d1_(1, 0) */
  id = map(1, 0, 1);
  d[1][id] = -d[1][map(1, 1, 0)];

  buildcoeff(ret[1], d[1], 1, a, r);

  for (i = 2; i <= l; i++) {
    for (m = 0; m <= i - 2; m++) {
      for (md = -m; md <= m; md++) {
        mid = map(i, m, md);

        c1 = (double)(i * (2 * i - 1)) /
             sqrt((double)((i * i - m * m) * (i * i - md * md)));

        c2 = d[1][4] - (m * md / (double)((i * (i - 1))));

        c3 = ((i - 1) * (i - 1) - m * m) * ((i - 1) * (i - 1) - md * md);
        c3 = sqrt(c3);
        c3 = c3 / (double)((i - 1) * (2 * i - 1));

        mid1 = map(i - 1, m, md);
        mid2 = map(i - 2, m, md);
        d[i][mid] = c1 * (c2 * d[i - 1][mid1] - c3 * d[i - 2][mid2]);
      }
    }

    /* dl_(l, l)		*/
    mid = map(i, i, i);
    mid1 = map(1, 1, 1);
    mid2 = map(i - 1, i - 1, i - 1);
    d[i][mid] = d[1][mid1] * d[i - 1][mid2];

    /* dl_(l-1, l-1)	*/
    mid = map(i, i - 1, i - 1);
    mid1 = map(1, 0, 0);
    mid2 = map(i - 1, i - 1, i - 1);

    d[i][mid] = (i * d[1][mid1] - i + 1) * d[i - 1][mid2];

    /* dl_(l, m - 1)	*/
    for (md = i - 1; md >= -i; md--) {
      mid = map(i, i, md);
      c1 = -sqrt((double)(i + (md + 1)) / (double)(i - (md + 1) + 1));

      d[i][mid] = c1 * tan(0.5 * b) * d[i][mid + 1];
    }

    /* dl_(l - 1, m - 1)	*/
    for (md = i - 2; md >= 1 - i; md--) {
      mid = map(i, i - 1, md);
      c1 = -(double)(i * cos(b) - (md + 1) + 1) /
           (double)(i * cos(b) - (md + 1));

      c1 *= sqrt((double)(i + (md + 1)) / (double)(i - (md + 1) + 1));

      d[i][mid] = c1 * tan(0.5 * b) * d[i][mid + 1];
    }

    fill_d(d[i], i);
    buildcoeff(ret[i], d[i], i, a, r);
  }

  for (i = 0; i <= l; i++) {
    free(d[i]);
  }
  free(d);
}

void buildcoeff(double *ret, double *c, int l, double a, double r) {
  int m, md;
  int id;
  int id1, id2;
  double w;

  for (m = -l; m <= l; m++) {
    for (md = -l; md <= l; md++) {
      id = map(l, m, md);
      id1 = (abs(md) + l) * (2 * l + 1) + (abs(m) + l);
      id2 = (abs(m) + l) * (2 * l + 1) + (-abs(md) + l);

      if (m == 0) {
        w = 1.0;
      } else {
        w = ((m % 2) == 0) ? 1.0 : -1.0;
      }

      ret[id] =
          sign(md) * phi(m, a) * phi(md, r) * ((c[id1] + w * c[id2]) * 0.5);

      ret[id] -=
          sign(m) * phi(-m, a) * phi(-md, r) * ((c[id1] - w * c[id2]) * 0.5);
    }
  }
}

void fill_d(double *m, int l) {
  int i, j;
  double w;

  /* left triangle of matrix with m>= 0 */
  for (j = 0; j < l; j++) {
    for (i = -l; i < -j; i++) {
      /* dl_(j, i) = dl_(-i, -j) */
      m[map(l, j, i)] = m[map(l, -i, -j)];
    }
  }

  /* right triangle of matrix with m>= 0 */
  for (j = 0; j < l; j++) {
    for (i = j + 1; i <= l; i++) {
      /* dl_(j, i) = (-1)^(j+i) dl_(i, j) */
      if (j + i == 0) {
        w = 1.0;
      } else {
        w = ((j + i) % 2 == 0) ? 1.0 : -1.0;
      }
      m[map(l, j, i)] = w * m[map(l, i, j)];
    }
  }
}

#endif /* TINYSHROT_IMPLEMENTATION_DEIFNED */
#endif /* TINYSHROT_IMPLEMENTATION */
