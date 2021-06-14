/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFT_for_cpp.h
 *
 * Code generation for function 'FFT_for_cpp'
 *
 */

#ifndef FFT_FOR_CPP_H
#define FFT_FOR_CPP_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "FFT_for_cpp_types.h"

/* Type Definitions */
class quick_fft
{
 public:
  quick_fft();
  ~quick_fft();
  static void FFT_for_cpp(const emxArray_real_T *inputData, real_T SampleFrequency,
                   emxArray_real_T *decibel, emxArray_real_T *Fv);

  static real_T FrequencyDetection(const emxArray_real_T *inputData, real_T
    SampleFrequency);
};

#endif

/* End of code generation (FFT_for_cpp.h) */
