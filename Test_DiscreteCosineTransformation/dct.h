//
//  dct.h
//  Test_DiscreteCosineTransformation
//
//  Created by Don on 10/21/15.
//  Copyright © 2015 Don. All rights reserved.
//

#ifndef dct_h
#define dct_h

#include <stdio.h>
#include <math.h>

void idct_ii_don(int N, const double X[], double x[]);

void dct_ii_don(unsigned char N,const double x[], double X[] );

// Implementation of LLM DCT.
void llm_dct(const double in[8], double out[8]);

void aan_dct(const double i[], double o[]);
#endif /* dct_h */

