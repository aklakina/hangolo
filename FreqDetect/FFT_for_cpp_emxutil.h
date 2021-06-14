/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFT_for_cpp_emxutil.h
 *
 * Code generation for function 'FFT_for_cpp_emxutil'
 *
 */

#ifndef FFT_FOR_CPP_EMXUTIL_H
#define FFT_FOR_CPP_EMXUTIL_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "FFT_for_cpp_types.h"

/* Function Declarations */
extern void emxEnsureCapacity_creal_T(emxArray_creal_T *emxArray, int32_T
  oldNumel);
extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int32_T
  oldNumel);
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int32_T oldNumel);
extern void emxFree_creal_T(emxArray_creal_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);

#endif

/* End of code generation (FFT_for_cpp_emxutil.h) */
