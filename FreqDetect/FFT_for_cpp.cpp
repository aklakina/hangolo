/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFT_for_cpp.cpp
 *
 * Code generation for function 'FFT_for_cpp'
 *
 */

/* Include files */
#include "FFT_for_cpp.h"
#include "FFT_for_cpp_emxutil.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <math.h>

/* Function Declarations */
static int32_T div_s32_floor(int32_T numerator, int32_T denominator);
static void dobluesteinfft(const emxArray_real_T *x, int32_T N2, int32_T n1,
  const emxArray_real_T *costab, const emxArray_real_T *sintab, const
  emxArray_real_T *sintabinv, emxArray_creal_T *y);
static void fft(const emxArray_real_T *x, real_T varargin_1, emxArray_creal_T *y);
static void r2br_r2dit_trig_impl(const emxArray_creal_T *x, int32_T
  unsigned_nRows, const emxArray_real_T *costab, const emxArray_real_T *sintab,
  emxArray_creal_T *y);
static real_T rt_hypotd_snf(real_T u0, real_T u1);
static real_T rt_powd_snf(real_T u0, real_T u1);

/* Function Definitions */
void quick_fft::FFT_for_cpp(const emxArray_real_T *inputData, real_T
  SampleFrequency, emxArray_real_T *decibel, emxArray_real_T *Fv)
{
  real_T n;
  int32_T eint;
  int32_T i;
  real_T delta1;
  emxArray_int32_T *Iv;
  emxArray_creal_T *x;
  int32_T nx;
  emxArray_real_T *b_decibel;
  n = frexp(static_cast<real_T>(inputData->size[1]), &eint);
  if (n == 0.5) {
    eint--;
  }

  n = rt_powd_snf(2.0, static_cast<real_T>(eint));

  /* point of FFT */
  i = Fv->size[0] * Fv->size[1];
  Fv->size[0] = 1;
  eint = static_cast<int32_T>(std::floor(n / 2.0));
  Fv->size[1] = eint;
  emxEnsureCapacity_real_T(Fv, i);
  if (eint >= 1) {
    Fv->data[eint - 1] = 1.0;
    if (Fv->size[1] >= 2) {
      Fv->data[0] = 0.0;
      if (Fv->size[1] >= 3) {
        delta1 = 1.0 / (static_cast<real_T>(Fv->size[1]) - 1.0);
        i = Fv->size[1];
        for (eint = 0; eint <= i - 3; eint++) {
          Fv->data[eint + 1] = (static_cast<real_T>(eint) + 1.0) * delta1;
        }
      }
    }
  }

  i = Fv->size[0] * Fv->size[1];
  eint = Fv->size[0] * Fv->size[1];
  Fv->size[0] = 1;
  emxEnsureCapacity_real_T(Fv, eint);
  eint = i - 1;
  for (i = 0; i <= eint; i++) {
    Fv->data[i] = Fv->data[i] * SampleFrequency / 2.0;
  }

  emxInit_int32_T(&Iv, 2);
  if (Fv->size[1] < 1) {
    Iv->size[0] = 1;
    Iv->size[1] = 0;
  } else {
    eint = Fv->size[1];
    i = Iv->size[0] * Iv->size[1];
    Iv->size[0] = 1;
    Iv->size[1] = Fv->size[1];
    emxEnsureCapacity_int32_T(Iv, i);
    for (i = 0; i < eint; i++) {
      Iv->data[i] = i + 1;
    }
  }

  emxInit_creal_T(&x, 2);
  fft(inputData, n, x);
  nx = x->size[1];
  i = decibel->size[0] * decibel->size[1];
  decibel->size[0] = 1;
  decibel->size[1] = x->size[1];
  emxEnsureCapacity_real_T(decibel, i);
  for (eint = 0; eint < nx; eint++) {
    decibel->data[eint] = rt_hypotd_snf(x->data[eint].re, x->data[eint].im);
  }

  emxFree_creal_T(&x);
  i = decibel->size[0] * decibel->size[1];
  eint = decibel->size[0] * decibel->size[1];
  decibel->size[0] = 1;
  emxEnsureCapacity_real_T(decibel, eint);
  eint = i - 1;
  for (i = 0; i <= eint; i++) {
    decibel->data[i] /= n;
  }

  emxInit_real_T(&b_decibel, 2);

  /*  calculate magnitude spectrum of Signal */
  i = b_decibel->size[0] * b_decibel->size[1];
  b_decibel->size[0] = 1;
  b_decibel->size[1] = Iv->size[1];
  emxEnsureCapacity_real_T(b_decibel, i);
  eint = Iv->size[0] * Iv->size[1];
  for (i = 0; i < eint; i++) {
    b_decibel->data[i] = decibel->data[Iv->data[i] - 1];
  }

  emxFree_int32_T(&Iv);
  i = decibel->size[0] * decibel->size[1];
  decibel->size[0] = 1;
  decibel->size[1] = b_decibel->size[1];
  emxEnsureCapacity_real_T(decibel, i);
  eint = b_decibel->size[0] * b_decibel->size[1];
  for (i = 0; i < eint; i++) {
    decibel->data[i] = b_decibel->data[i];
  }

  emxFree_real_T(&b_decibel);

  /*  One-sided spectrum */
  nx = decibel->size[1];
  for (eint = 0; eint < nx; eint++) {
    decibel->data[eint] = std::log10(decibel->data[eint]);
  }

  i = decibel->size[0] * decibel->size[1];
  eint = decibel->size[0] * decibel->size[1];
  decibel->size[0] = 1;
  emxEnsureCapacity_real_T(decibel, eint);
  eint = i - 1;
  for (i = 0; i <= eint; i++) {
    decibel->data[i] *= 20.0;
  }

  /* Convert to decibels */
}

real_T quick_fft::FrequencyDetection(const emxArray_real_T *inputData, real_T
  SampleFrequency)
{
  real_T fundFreq;
  real_T n;
  int32_T nx;
  emxArray_real_T *Fv;
  int32_T i;
  int32_T k;
  real_T delta1;
  emxArray_real_T *Iv;
  emxArray_real_T *b_y0;
  emxArray_creal_T *x;
  emxArray_real_T *c_y0;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  int32_T i6;
  int32_T sizeX_idx_1_tmp;
  int32_T b_n;
  boolean_T exitg1;
  n = frexp(static_cast<real_T>(inputData->size[1]), &nx);
  if (n == 0.5) {
    nx--;
  }

  emxInit_real_T(&Fv, 2);
  n = rt_powd_snf(2.0, static_cast<real_T>(nx));

  /* point of FFT */
  i = Fv->size[0] * Fv->size[1];
  Fv->size[0] = 1;
  k = static_cast<int32_T>(std::floor(n / 2.0));
  Fv->size[1] = k;
  emxEnsureCapacity_real_T(Fv, i);
  if (k >= 1) {
    Fv->data[k - 1] = 1.0;
    if (Fv->size[1] >= 2) {
      Fv->data[0] = 0.0;
      if (Fv->size[1] >= 3) {
        delta1 = 1.0 / (static_cast<real_T>(Fv->size[1]) - 1.0);
        i = Fv->size[1];
        for (k = 0; k <= i - 3; k++) {
          Fv->data[k + 1] = (static_cast<real_T>(k) + 1.0) * delta1;
        }
      }
    }
  }

  i = Fv->size[0] * Fv->size[1];
  k = Fv->size[0] * Fv->size[1];
  Fv->size[0] = 1;
  emxEnsureCapacity_real_T(Fv, k);
  nx = i - 1;
  for (i = 0; i <= nx; i++) {
    Fv->data[i] = Fv->data[i] * SampleFrequency / 2.0;
  }

  emxInit_real_T(&Iv, 2);
  if (Fv->size[1] < 1) {
    Iv->size[0] = 1;
    Iv->size[1] = 0;
  } else {
    i = Iv->size[0] * Iv->size[1];
    Iv->size[0] = 1;
    Iv->size[1] = static_cast<int32_T>(static_cast<real_T>(Fv->size[1]) - 1.0) +
      1;
    emxEnsureCapacity_real_T(Iv, i);
    nx = static_cast<int32_T>(static_cast<real_T>(Fv->size[1]) - 1.0);
    for (i = 0; i <= nx; i++) {
      Iv->data[i] = static_cast<real_T>(i) + 1.0;
    }
  }

  emxInit_real_T(&b_y0, 2);
  emxInit_creal_T(&x, 2);
  fft(inputData, n, x);
  nx = x->size[1];
  i = b_y0->size[0] * b_y0->size[1];
  b_y0->size[0] = 1;
  b_y0->size[1] = x->size[1];
  emxEnsureCapacity_real_T(b_y0, i);
  for (k = 0; k < nx; k++) {
    b_y0->data[k] = rt_hypotd_snf(x->data[k].re, x->data[k].im);
  }

  emxFree_creal_T(&x);
  i = b_y0->size[0] * b_y0->size[1];
  k = b_y0->size[0] * b_y0->size[1];
  b_y0->size[0] = 1;
  emxEnsureCapacity_real_T(b_y0, k);
  nx = i - 1;
  for (i = 0; i <= nx; i++) {
    b_y0->data[i] /= n;
  }

  emxInit_real_T(&c_y0, 2);

  /*  calculate magnitude spectrum of Signal */
  i = c_y0->size[0] * c_y0->size[1];
  c_y0->size[0] = 1;
  c_y0->size[1] = Iv->size[1];
  emxEnsureCapacity_real_T(c_y0, i);
  nx = Iv->size[0] * Iv->size[1];
  for (i = 0; i < nx; i++) {
    c_y0->data[i] = b_y0->data[static_cast<int32_T>(Iv->data[i]) - 1];
  }

  i = b_y0->size[0] * b_y0->size[1];
  b_y0->size[0] = 1;
  b_y0->size[1] = c_y0->size[1];
  emxEnsureCapacity_real_T(b_y0, i);
  nx = c_y0->size[0] * c_y0->size[1];
  for (i = 0; i < nx; i++) {
    b_y0->data[i] = c_y0->data[i];
  }

  emxFree_real_T(&c_y0);

  /*  One-sided spectrum */
  /* Harmonic Product Spectrum */
  if (1 > b_y0->size[1]) {
    i = 0;
  } else {
    i = b_y0->size[1];
  }

  if (1 > b_y0->size[1]) {
    k = 1;
    i1 = 0;
  } else {
    k = 2;
    i1 = b_y0->size[1];
  }

  if (1 > b_y0->size[1]) {
    i2 = 1;
    i3 = 0;
  } else {
    i2 = 3;
    i3 = b_y0->size[1];
  }

  if (1 > b_y0->size[1]) {
    i4 = 1;
    i5 = 0;
  } else {
    i4 = 4;
    i5 = b_y0->size[1];
  }

  if (1 > b_y0->size[1]) {
    i6 = 1;
    nx = -1;
  } else {
    i6 = 5;
    nx = b_y0->size[1] - 1;
  }

  sizeX_idx_1_tmp = div_s32_floor(nx, i6) + 1;
  if (sizeX_idx_1_tmp == 0) {
    b_n = 0;
  } else {
    b_n = sizeX_idx_1_tmp;
  }

  nx = Iv->size[0] * Iv->size[1];
  Iv->size[0] = 1;
  Iv->size[1] = b_n;
  emxEnsureCapacity_real_T(Iv, nx);
  for (nx = 0; nx < b_n; nx++) {
    Iv->data[nx] = 0.0;
  }

  if (sizeX_idx_1_tmp == 0) {
    b_n = 0;
  } else {
    b_n = sizeX_idx_1_tmp;
  }

  for (nx = 0; nx < b_n; nx++) {
    Iv->data[nx] = b_y0->data[nx % i] * b_y0->data[k * (nx % (div_s32_floor(i1 -
      1, k) + 1))] * b_y0->data[i2 * (nx % (div_s32_floor(i3 - 1, i2) + 1))] *
      b_y0->data[i4 * (nx % (div_s32_floor(i5 - 1, i4) + 1))] * b_y0->data[i6 *
      (nx % sizeX_idx_1_tmp)];
  }

  emxFree_real_T(&b_y0);

  /* Finding the fundamental frequency */
  b_n = Iv->size[1];
  if (Iv->size[1] <= 2) {
    if (Iv->size[1] == 1) {
      nx = 1;
    } else if ((Iv->data[0] < Iv->data[1]) || (rtIsNaN(Iv->data[0]) && (!rtIsNaN
                 (Iv->data[1])))) {
      nx = 2;
    } else {
      nx = 1;
    }
  } else {
    if (!rtIsNaN(Iv->data[0])) {
      nx = 1;
    } else {
      nx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= Iv->size[1])) {
        if (!rtIsNaN(Iv->data[k - 1])) {
          nx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (nx == 0) {
      nx = 1;
    } else {
      n = Iv->data[nx - 1];
      i = nx + 1;
      for (k = i; k <= b_n; k++) {
        delta1 = Iv->data[k - 1];
        if (n < delta1) {
          n = delta1;
          nx = k;
        }
      }
    }
  }

  emxFree_real_T(&Iv);
  fundFreq = Fv->data[nx - 1];
  emxFree_real_T(&Fv);
  return fundFreq;
}

quick_fft::~quick_fft()
{
  /* (no terminate code required) */
}

quick_fft::quick_fft()
{
  rt_InitInfAndNaN();
}

static int32_T div_s32_floor(int32_T numerator, int32_T denominator)
{
  int32_T quotient;
  uint32_T absNumerator;
  uint32_T absDenominator;
  boolean_T quotientNeedsNegation;
  uint32_T tempAbsQuotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~static_cast<uint32_T>(numerator) + 1U;
    } else {
      absNumerator = static_cast<uint32_T>(numerator);
    }

    if (denominator < 0) {
      absDenominator = ~static_cast<uint32_T>(denominator) + 1U;
    } else {
      absDenominator = static_cast<uint32_T>(denominator);
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if (quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > 0U) {
        tempAbsQuotient++;
      }

      quotient = -static_cast<int32_T>(tempAbsQuotient);
    } else {
      quotient = static_cast<int32_T>(tempAbsQuotient);
    }
  }

  return quotient;
}

static void dobluesteinfft(const emxArray_real_T *x, int32_T N2, int32_T n1,
  const emxArray_real_T *costab, const emxArray_real_T *sintab, const
  emxArray_real_T *sintabinv, emxArray_creal_T *y)
{
  emxArray_creal_T *wwc;
  int32_T nInt2m1;
  int32_T ju;
  int32_T idx;
  int32_T rt;
  int32_T nInt2;
  int32_T ihi;
  real_T nt_im;
  real_T nt_re;
  int32_T nt_re_tmp;
  emxArray_creal_T *fy;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T i;
  boolean_T tst;
  real_T twid_re;
  real_T twid_im;
  emxArray_creal_T *fv;
  emxInit_creal_T(&wwc, 1);
  nInt2m1 = (n1 + n1) - 1;
  ju = wwc->size[0];
  wwc->size[0] = nInt2m1;
  emxEnsureCapacity_creal_T(wwc, ju);
  idx = n1;
  rt = 0;
  wwc->data[n1 - 1].re = 1.0;
  wwc->data[n1 - 1].im = 0.0;
  nInt2 = n1 << 1;
  for (ihi = 0; ihi <= n1 - 2; ihi++) {
    ju = ((ihi + 1) << 1) - 1;
    if (nInt2 - rt <= ju) {
      rt += ju - nInt2;
    } else {
      rt += ju;
    }

    nt_im = -3.1415926535897931 * static_cast<real_T>(rt) / static_cast<real_T>
      (n1);
    if (nt_im == 0.0) {
      nt_re = 1.0;
      nt_im = 0.0;
    } else {
      nt_re = std::cos(nt_im);
      nt_im = std::sin(nt_im);
    }

    wwc->data[idx - 2].re = nt_re;
    wwc->data[idx - 2].im = -nt_im;
    idx--;
  }

  idx = 0;
  ju = nInt2m1 - 1;
  for (ihi = ju; ihi >= n1; ihi--) {
    wwc->data[ihi] = wwc->data[idx];
    idx++;
  }

  nInt2 = x->size[0];
  if (n1 < nInt2) {
    nInt2 = n1;
  }

  ju = y->size[0];
  y->size[0] = n1;
  emxEnsureCapacity_creal_T(y, ju);
  if (n1 > x->size[0]) {
    ju = y->size[0];
    y->size[0] = n1;
    emxEnsureCapacity_creal_T(y, ju);
    for (ju = 0; ju < n1; ju++) {
      y->data[ju].re = 0.0;
      y->data[ju].im = 0.0;
    }
  }

  rt = 0;
  for (ihi = 0; ihi < nInt2; ihi++) {
    nt_re_tmp = (n1 + ihi) - 1;
    y->data[ihi].re = wwc->data[nt_re_tmp].re * x->data[rt];
    y->data[ihi].im = wwc->data[nt_re_tmp].im * -x->data[rt];
    rt++;
  }

  ju = nInt2 + 1;
  for (ihi = ju; ihi <= n1; ihi++) {
    y->data[ihi - 1].re = 0.0;
    y->data[ihi - 1].im = 0.0;
  }

  emxInit_creal_T(&fy, 1);
  nInt2m1 = y->size[0];
  if (nInt2m1 >= N2) {
    nInt2m1 = N2;
  }

  idx = N2 - 2;
  nRowsD2 = N2 / 2;
  nRowsD4 = nRowsD2 / 2;
  ju = fy->size[0];
  fy->size[0] = N2;
  emxEnsureCapacity_creal_T(fy, ju);
  if (N2 > y->size[0]) {
    ju = fy->size[0];
    fy->size[0] = N2;
    emxEnsureCapacity_creal_T(fy, ju);
    for (ju = 0; ju < N2; ju++) {
      fy->data[ju].re = 0.0;
      fy->data[ju].im = 0.0;
    }
  }

  nInt2 = 0;
  ju = 0;
  rt = 0;
  for (i = 0; i <= nInt2m1 - 2; i++) {
    fy->data[rt] = y->data[nInt2];
    rt = N2;
    tst = true;
    while (tst) {
      rt >>= 1;
      ju ^= rt;
      tst = ((ju & rt) == 0);
    }

    rt = ju;
    nInt2++;
  }

  fy->data[rt] = y->data[nInt2];
  if (N2 > 1) {
    for (i = 0; i <= idx; i += 2) {
      nt_re = fy->data[i + 1].re;
      nt_im = fy->data[i + 1].im;
      twid_re = fy->data[i].re;
      twid_im = fy->data[i].im;
      fy->data[i + 1].re = fy->data[i].re - fy->data[i + 1].re;
      fy->data[i + 1].im = fy->data[i].im - fy->data[i + 1].im;
      twid_re += nt_re;
      twid_im += nt_im;
      fy->data[i].re = twid_re;
      fy->data[i].im = twid_im;
    }
  }

  rt = 2;
  nInt2 = 4;
  nInt2m1 = ((nRowsD4 - 1) << 2) + 1;
  while (nRowsD4 > 0) {
    for (i = 0; i < nInt2m1; i += nInt2) {
      nt_re_tmp = i + rt;
      nt_re = fy->data[nt_re_tmp].re;
      nt_im = fy->data[nt_re_tmp].im;
      fy->data[nt_re_tmp].re = fy->data[i].re - fy->data[nt_re_tmp].re;
      fy->data[nt_re_tmp].im = fy->data[i].im - fy->data[nt_re_tmp].im;
      fy->data[i].re += nt_re;
      fy->data[i].im += nt_im;
    }

    idx = 1;
    for (ju = nRowsD4; ju < nRowsD2; ju += nRowsD4) {
      twid_re = costab->data[ju];
      twid_im = sintab->data[ju];
      i = idx;
      ihi = idx + nInt2m1;
      while (i < ihi) {
        nt_re_tmp = i + rt;
        nt_re = twid_re * fy->data[nt_re_tmp].re - twid_im * fy->data[nt_re_tmp]
          .im;
        nt_im = twid_re * fy->data[nt_re_tmp].im + twid_im * fy->data[nt_re_tmp]
          .re;
        fy->data[nt_re_tmp].re = fy->data[i].re - nt_re;
        fy->data[nt_re_tmp].im = fy->data[i].im - nt_im;
        fy->data[i].re += nt_re;
        fy->data[i].im += nt_im;
        i += nInt2;
      }

      idx++;
    }

    nRowsD4 /= 2;
    rt = nInt2;
    nInt2 += nInt2;
    nInt2m1 -= rt;
  }

  emxInit_creal_T(&fv, 1);
  r2br_r2dit_trig_impl(wwc, N2, costab, sintab, fv);
  rt = fy->size[0];
  for (ju = 0; ju < rt; ju++) {
    twid_im = fy->data[ju].re * fv->data[ju].im + fy->data[ju].im * fv->data[ju]
      .re;
    fy->data[ju].re = fy->data[ju].re * fv->data[ju].re - fy->data[ju].im *
      fv->data[ju].im;
    fy->data[ju].im = twid_im;
  }

  r2br_r2dit_trig_impl(fy, N2, costab, sintabinv, fv);
  emxFree_creal_T(&fy);
  if (fv->size[0] > 1) {
    twid_re = 1.0 / static_cast<real_T>(fv->size[0]);
    rt = fv->size[0];
    for (ju = 0; ju < rt; ju++) {
      fv->data[ju].re *= twid_re;
      fv->data[ju].im *= twid_re;
    }
  }

  idx = 0;
  ju = wwc->size[0];
  for (ihi = n1; ihi <= ju; ihi++) {
    y->data[idx].re = wwc->data[ihi - 1].re * fv->data[ihi - 1].re + wwc->
      data[ihi - 1].im * fv->data[ihi - 1].im;
    y->data[idx].im = wwc->data[ihi - 1].re * fv->data[ihi - 1].im - wwc->
      data[ihi - 1].im * fv->data[ihi - 1].re;
    idx++;
  }

  emxFree_creal_T(&fv);
  emxFree_creal_T(&wwc);
}

static void fft(const emxArray_real_T *x, real_T varargin_1, emxArray_creal_T *y)
{
  emxArray_creal_T *b_y;
  emxArray_real_T *costab1q;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  emxArray_real_T *sintabinv;
  boolean_T guard1 = false;
  int32_T loop_ub;
  boolean_T useRadix2;
  int32_T pmin;
  int32_T pow2p;
  int32_T ihi;
  int32_T pmax;
  real_T twid_re;
  boolean_T exitg1;
  int32_T istart;
  emxArray_real_T b_x;
  int32_T c_x[1];
  int32_T nRowsM2;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T i;
  real_T temp_re;
  real_T temp_im;
  real_T twid_im;
  emxInit_creal_T(&b_y, 1);
  emxInit_real_T(&costab1q, 2);
  emxInit_real_T(&costab, 2);
  emxInit_real_T(&sintab, 2);
  emxInit_real_T(&sintabinv, 2);
  guard1 = false;
  if (x->size[1] == 0) {
    guard1 = true;
  } else {
    loop_ub = static_cast<int32_T>(varargin_1);
    if (loop_ub == 0) {
      guard1 = true;
    } else {
      useRadix2 = ((loop_ub & (loop_ub - 1)) == 0);
      pmin = 1;
      if (useRadix2) {
        pmax = loop_ub;
      } else {
        ihi = (loop_ub + loop_ub) - 1;
        pmax = 31;
        if (ihi <= 1) {
          pmax = 0;
        } else {
          pmin = 0;
          exitg1 = false;
          while ((!exitg1) && (pmax - pmin > 1)) {
            istart = (pmin + pmax) >> 1;
            pow2p = 1 << istart;
            if (pow2p == ihi) {
              pmax = istart;
              exitg1 = true;
            } else if (pow2p > ihi) {
              pmax = istart;
            } else {
              pmin = istart;
            }
          }
        }

        pmin = 1 << pmax;
        pmax = pmin;
      }

      twid_re = 6.2831853071795862 / static_cast<real_T>(pmax);
      ihi = pmax / 2 / 2;
      pow2p = costab1q->size[0] * costab1q->size[1];
      costab1q->size[0] = 1;
      costab1q->size[1] = ihi + 1;
      emxEnsureCapacity_real_T(costab1q, pow2p);
      costab1q->data[0] = 1.0;
      pmax = ihi / 2 - 1;
      for (istart = 0; istart <= pmax; istart++) {
        costab1q->data[istart + 1] = std::cos(twid_re * (static_cast<real_T>
          (istart) + 1.0));
      }

      pow2p = pmax + 2;
      pmax = ihi - 1;
      for (istart = pow2p; istart <= pmax; istart++) {
        costab1q->data[istart] = std::sin(twid_re * static_cast<real_T>(ihi -
          istart));
      }

      costab1q->data[ihi] = 0.0;
      if (!useRadix2) {
        ihi = costab1q->size[1] - 1;
        pmax = (costab1q->size[1] - 1) << 1;
        pow2p = costab->size[0] * costab->size[1];
        costab->size[0] = 1;
        costab->size[1] = pmax + 1;
        emxEnsureCapacity_real_T(costab, pow2p);
        pow2p = sintab->size[0] * sintab->size[1];
        sintab->size[0] = 1;
        sintab->size[1] = pmax + 1;
        emxEnsureCapacity_real_T(sintab, pow2p);
        costab->data[0] = 1.0;
        sintab->data[0] = 0.0;
        pow2p = sintabinv->size[0] * sintabinv->size[1];
        sintabinv->size[0] = 1;
        sintabinv->size[1] = pmax + 1;
        emxEnsureCapacity_real_T(sintabinv, pow2p);
        for (istart = 0; istart < ihi; istart++) {
          sintabinv->data[istart + 1] = costab1q->data[(ihi - istart) - 1];
        }

        pow2p = costab1q->size[1];
        for (istart = pow2p; istart <= pmax; istart++) {
          sintabinv->data[istart] = costab1q->data[istart - ihi];
        }

        for (istart = 0; istart < ihi; istart++) {
          costab->data[istart + 1] = costab1q->data[istart + 1];
          sintab->data[istart + 1] = -costab1q->data[(ihi - istart) - 1];
        }

        pow2p = costab1q->size[1];
        for (istart = pow2p; istart <= pmax; istart++) {
          costab->data[istart] = -costab1q->data[pmax - istart];
          sintab->data[istart] = -costab1q->data[istart - ihi];
        }
      } else {
        ihi = costab1q->size[1] - 1;
        pmax = (costab1q->size[1] - 1) << 1;
        pow2p = costab->size[0] * costab->size[1];
        costab->size[0] = 1;
        costab->size[1] = pmax + 1;
        emxEnsureCapacity_real_T(costab, pow2p);
        pow2p = sintab->size[0] * sintab->size[1];
        sintab->size[0] = 1;
        sintab->size[1] = pmax + 1;
        emxEnsureCapacity_real_T(sintab, pow2p);
        costab->data[0] = 1.0;
        sintab->data[0] = 0.0;
        for (istart = 0; istart < ihi; istart++) {
          costab->data[istart + 1] = costab1q->data[istart + 1];
          sintab->data[istart + 1] = -costab1q->data[(ihi - istart) - 1];
        }

        pow2p = costab1q->size[1];
        for (istart = pow2p; istart <= pmax; istart++) {
          costab->data[istart] = -costab1q->data[pmax - istart];
          sintab->data[istart] = -costab1q->data[istart - ihi];
        }

        sintabinv->size[0] = 1;
        sintabinv->size[1] = 0;
      }

      if (useRadix2) {
        istart = x->size[1];
        if (istart >= loop_ub) {
          istart = loop_ub;
        }

        nRowsM2 = loop_ub - 2;
        nRowsD2 = loop_ub / 2;
        nRowsD4 = nRowsD2 / 2;
        pow2p = b_y->size[0];
        b_y->size[0] = loop_ub;
        emxEnsureCapacity_creal_T(b_y, pow2p);
        if (loop_ub > x->size[1]) {
          pow2p = b_y->size[0];
          b_y->size[0] = loop_ub;
          emxEnsureCapacity_creal_T(b_y, pow2p);
          for (pow2p = 0; pow2p < loop_ub; pow2p++) {
            b_y->data[pow2p].re = 0.0;
            b_y->data[pow2p].im = 0.0;
          }
        }

        pmin = 0;
        pow2p = 0;
        pmax = 0;
        for (i = 0; i <= istart - 2; i++) {
          b_y->data[pmax].re = x->data[pmin];
          b_y->data[pmax].im = 0.0;
          ihi = loop_ub;
          useRadix2 = true;
          while (useRadix2) {
            ihi >>= 1;
            pow2p ^= ihi;
            useRadix2 = ((pow2p & ihi) == 0);
          }

          pmax = pow2p;
          pmin++;
        }

        b_y->data[pmax].re = x->data[pmin];
        b_y->data[pmax].im = 0.0;
        if (loop_ub > 1) {
          for (i = 0; i <= nRowsM2; i += 2) {
            temp_re = b_y->data[i + 1].re;
            temp_im = b_y->data[i + 1].im;
            twid_re = b_y->data[i].re;
            twid_im = b_y->data[i].im;
            b_y->data[i + 1].re = b_y->data[i].re - b_y->data[i + 1].re;
            b_y->data[i + 1].im = b_y->data[i].im - b_y->data[i + 1].im;
            twid_re += temp_re;
            twid_im += temp_im;
            b_y->data[i].re = twid_re;
            b_y->data[i].im = twid_im;
          }
        }

        pmax = 2;
        pmin = 4;
        pow2p = ((nRowsD4 - 1) << 2) + 1;
        while (nRowsD4 > 0) {
          for (i = 0; i < pow2p; i += pmin) {
            loop_ub = i + pmax;
            temp_re = b_y->data[loop_ub].re;
            temp_im = b_y->data[loop_ub].im;
            b_y->data[loop_ub].re = b_y->data[i].re - b_y->data[loop_ub].re;
            b_y->data[loop_ub].im = b_y->data[i].im - b_y->data[loop_ub].im;
            b_y->data[i].re += temp_re;
            b_y->data[i].im += temp_im;
          }

          istart = 1;
          for (nRowsM2 = nRowsD4; nRowsM2 < nRowsD2; nRowsM2 += nRowsD4) {
            twid_re = costab->data[nRowsM2];
            twid_im = sintab->data[nRowsM2];
            i = istart;
            ihi = istart + pow2p;
            while (i < ihi) {
              loop_ub = i + pmax;
              temp_re = twid_re * b_y->data[loop_ub].re - twid_im * b_y->
                data[loop_ub].im;
              temp_im = twid_re * b_y->data[loop_ub].im + twid_im * b_y->
                data[loop_ub].re;
              b_y->data[loop_ub].re = b_y->data[i].re - temp_re;
              b_y->data[loop_ub].im = b_y->data[i].im - temp_im;
              b_y->data[i].re += temp_re;
              b_y->data[i].im += temp_im;
              i += pmin;
            }

            istart++;
          }

          nRowsD4 /= 2;
          pmax = pmin;
          pmin += pmin;
          pow2p -= pmax;
        }
      } else {
        pmax = x->size[1];
        b_x = *x;
        c_x[0] = pmax;
        b_x.size = &c_x[0];
        b_x.numDimensions = 1;
        dobluesteinfft(&b_x, pmin, loop_ub, costab, sintab, sintabinv, b_y);
      }
    }
  }

  if (guard1) {
    loop_ub = static_cast<int32_T>(varargin_1);
    pow2p = b_y->size[0];
    b_y->size[0] = loop_ub;
    emxEnsureCapacity_creal_T(b_y, pow2p);
    if (loop_ub > x->size[1]) {
      pow2p = b_y->size[0];
      b_y->size[0] = loop_ub;
      emxEnsureCapacity_creal_T(b_y, pow2p);
      for (pow2p = 0; pow2p < loop_ub; pow2p++) {
        b_y->data[pow2p].re = 0.0;
        b_y->data[pow2p].im = 0.0;
      }
    }
  }

  emxFree_real_T(&sintabinv);
  emxFree_real_T(&sintab);
  emxFree_real_T(&costab);
  emxFree_real_T(&costab1q);
  pow2p = y->size[0] * y->size[1];
  y->size[0] = 1;
  loop_ub = static_cast<int32_T>(varargin_1);
  y->size[1] = loop_ub;
  emxEnsureCapacity_creal_T(y, pow2p);
  for (pow2p = 0; pow2p < loop_ub; pow2p++) {
    y->data[pow2p] = b_y->data[pow2p];
  }

  emxFree_creal_T(&b_y);
}

static void r2br_r2dit_trig_impl(const emxArray_creal_T *x, int32_T
  unsigned_nRows, const emxArray_real_T *costab, const emxArray_real_T *sintab,
  emxArray_creal_T *y)
{
  int32_T istart;
  int32_T nRowsM2;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T iy;
  int32_T ix;
  int32_T ju;
  int32_T i;
  boolean_T tst;
  real_T temp_re;
  real_T temp_im;
  real_T twid_re;
  real_T twid_im;
  int32_T temp_re_tmp;
  int32_T ihi;
  istart = x->size[0];
  if (istart >= unsigned_nRows) {
    istart = unsigned_nRows;
  }

  nRowsM2 = unsigned_nRows - 2;
  nRowsD2 = unsigned_nRows / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity_creal_T(y, iy);
  if (unsigned_nRows > x->size[0]) {
    iy = y->size[0];
    y->size[0] = unsigned_nRows;
    emxEnsureCapacity_creal_T(y, iy);
    for (iy = 0; iy < unsigned_nRows; iy++) {
      y->data[iy].re = 0.0;
      y->data[iy].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i <= istart - 2; i++) {
    y->data[iy] = x->data[ix];
    iy = unsigned_nRows;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy] = x->data[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= nRowsM2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      twid_re = y->data[i].re;
      twid_im = y->data[i].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      twid_re += temp_re;
      twid_im += temp_im;
      y->data[i].re = twid_re;
      y->data[i].im = twid_im;
    }
  }

  iy = 2;
  ix = 4;
  ju = ((nRowsD4 - 1) << 2) + 1;
  while (nRowsD4 > 0) {
    for (i = 0; i < ju; i += ix) {
      temp_re_tmp = i + iy;
      temp_re = y->data[temp_re_tmp].re;
      temp_im = y->data[temp_re_tmp].im;
      y->data[temp_re_tmp].re = y->data[i].re - y->data[temp_re_tmp].re;
      y->data[temp_re_tmp].im = y->data[i].im - y->data[temp_re_tmp].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    istart = 1;
    for (nRowsM2 = nRowsD4; nRowsM2 < nRowsD2; nRowsM2 += nRowsD4) {
      twid_re = costab->data[nRowsM2];
      twid_im = sintab->data[nRowsM2];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re_tmp = i + iy;
        temp_re = twid_re * y->data[temp_re_tmp].re - twid_im * y->
          data[temp_re_tmp].im;
        temp_im = twid_re * y->data[temp_re_tmp].im + twid_im * y->
          data[temp_re_tmp].re;
        y->data[temp_re_tmp].re = y->data[i].re - temp_re;
        y->data[temp_re_tmp].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    nRowsD4 /= 2;
    iy = ix;
    ix += ix;
    ju -= iy;
  }
}

static real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = std::abs(u0);
  y = std::abs(u1);
  if (a < y) {
    a /= y;
    y *= std::sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * std::sqrt(y * y + 1.0);
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

static real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T d;
  real_T d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* End of code generation (FFT_for_cpp.cpp) */
