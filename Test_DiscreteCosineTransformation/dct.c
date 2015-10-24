//
//  dct.c
//  Test_DiscreteCosineTransformation
//
//  Created by Don on 10/21/15.
//  Copyright © 2015 Don. All rights reserved.
//

#include "dct.h"

void idct_ii_don(int N, const double X[], double x[]) {
    for (int n=0; n < N; ++n) {
        double sum = 0;
        for (int k=0; k<N; ++k) {
            double s = (k == 0) ? sqrt(.5) : 1.;
            sum += s * X[k] * cos(M_PI * (n + 0.5) * k / N);
        }
        x[n] = sum * sqrt(2. / N);
    }
}

void dct_ii_don(unsigned char N,const double x[], double X[] )
{
    for (unsigned char k = 0; k < N; ++k ) {
        double sum = 0.;
        double s = (k==0) ? sqrt(.5) : 1.;
        for (unsigned char n = 0; n < N; ++n) {
            sum += s * x[n] * cos(M_PI * (n+ 0.5) * k /N );
        }
        X[k] = sum * sqrt(2.0 / N);
    }
    //stp x29, x30, [sp, #-16]!
    //mov x29, sp
    //sub sp, sp, #80
    //sturb w0, [x29, #-1]   ;1st arg, unsigned char N
    //stur x1, [x29, #-16]  ;2nd arg, const double x[]
    //stur x2, [x29, #-24]  ;3rd arg, double X[]
    //sturb wzr, [x29, #-25];
    //ldurb w8, [x29, #-25]  ;unsigned char k = 0
    //ldurb w9, [x29, #-1]  ; N
    //cmp w8, w9
    //b.ge 0x....
    
    
}

#define SQRT_2 1.4142135623730950488016887242097
#define SQRT_8 2.8284271247461900976033774484194

// Implementation of LLM DCT.
void llm_dct(const double in[8], double out[8]) {
    // Constants:
//    const double s1 = sin(1. * M_PI / 16.);
//    const double c1 = cos(1. * M_PI / 16.);
//    const double s3 = sin(3. * M_PI / 16.);
//    const double c3 = cos(3. * M_PI / 16.);
//    const double r2s6 = sqrt(2.) * sin(6. * M_PI / 16.);
//    const double r2c6 = sqrt(2.) * cos(6. * M_PI / 16.);

    const double s1 = 0.19509032201612826784828486847702;
    const double c1 = 0.98078528040323044912618223613424;
    const double s3 = 0.55557023301960222474283081394853;
    const double c3 = 0.83146961230254523707878837761791;
    const double r2s6 = 1.4142135623730950488016887242097 * 0.92387953251128675612818318939679;
    const double r2c6 = 1.4142135623730950488016887242097 * 0.3826834323650897717284599840304;
    
    // After stage 1:
    const double s1_0 =  in[0] + in[7];
    const double s1_1 =  in[1] + in[6];
    const double s1_2 =  in[2] + in[5];
    const double s1_3 =  in[3] + in[4];
    const double s1_4 = -in[4] + in[3];
    const double s1_5 = -in[5] + in[2];
    const double s1_6 = -in[6] + in[1];
    const double s1_7 = -in[7] + in[0];
    
    // After stage 2:
    const double s2_0 =  s1_0 + s1_3;
    const double s2_1 =  s1_1 + s1_2;
    const double s2_2 = -s1_2 + s1_1;
    const double s2_3 = -s1_3 + s1_0;
    
    const double z1 = c3 * (s1_7 + s1_4);
    const double s2_4 = ( s3-c3) * s1_7 + z1;
    const double s2_7 = (-s3-c3) * s1_4 + z1;
    
    const double z2 = c1 * (s1_6 + s1_5);
    const double s2_5 = ( s1-c1) * s1_6 + z2;
    const double s2_6 = (-s1-c1) * s1_5 + z2;
    
    // After stage 3:
    const double s3_0 =  s2_0 + s2_1;
    const double s3_1 = -s2_1 + s2_0;
    
    const double z3 = r2c6 * (s2_3 + s2_2);
    const double s3_2 = ( r2s6-r2c6) * s2_3 + z3;
    const double s3_3 = (-r2s6-r2c6) * s2_2 + z3;
    
    const double s3_4 =  s2_4 + s2_6;
    const double s3_5 = -s2_5 + s2_7;
    const double s3_6 = -s2_6 + s2_4;
    const double s3_7 =  s2_7 + s2_5;
    
    // After stage 4:
    const double s4_4 = -s3_4 + s3_7;
    const double s4_5 =  s3_5 * SQRT_2;
    const double s4_6 =  s3_6 * SQRT_2;
    const double s4_7 =  s3_7 + s3_4;
    
    // Shuffle and scaling:
    out[0] = s3_0 / SQRT_8;
    out[4] = s3_1 / SQRT_8;
    out[2] = s3_2 / SQRT_8;
    out[6] = s3_3 / SQRT_8;
    out[7] = s4_4 / SQRT_8;
    out[3] = s4_5 / SQRT_8;  // Alternative: s3_5 / 2
    out[5] = s4_6 / SQRT_8;
    out[1] = s4_7 / SQRT_8;
}



void aan_dct(const double i[], double o[]) {
#if 1
    const double a1 = sqrt(.5);
    const double a2 = sqrt(2.) * cos(3./16. * 2 * M_PI);
    const double a3 = a1;
    const double a4 = sqrt(2.) * cos(1./16. * 2 * M_PI);
    const double a5 = cos(3./16. * 2 * M_PI);
#else
    const double a1 = 0.707;
    const double a2 = 0.541;
    const double a3 = 0.707;
    const double a4 = 1.307;
    const double a5 = 0.383;
#endif
    
    double b0 = i[0] + i[7];
    double b1 = i[1] + i[6];
    double b2 = i[2] + i[5];
    double b3 = i[3] + i[4];
    double b4 =-i[4] + i[3];
    double b5 =-i[5] + i[2];
    double b6 =-i[6] + i[1];
    double b7 =-i[7] + i[0];
    
    double c0 = b0 + b3;
    double c1 = b1 + b2;
    double c2 =-b2 + b1;
    double c3 =-b3 + b0;
    double c4 =-b4 - b5;
    double c5 = b5 + b6;
    double c6 = b6 + b7;
    double c7 = b7;
    
    double d0 = c0 + c1;
    double d1 =-c1 + c0;
    double d2 = c2 + c3;
    double d3 = c3;
    double d4 = c4;
    double d5 = c5;
    double d6 = c6;
    double d7 = c7;
    
    double d8 = (d4 + d6) * a5;
    
    double e0 = d0;
    double e1 = d1;
    double e2 = d2 * a1;
    double e3 = d3;
    double e4 = -d4 * a2 - d8;
    double e5 = d5 * a3;
    double e6 = d6 * a4 - d8;
    double e7 = d7;
    
    double f0 = e0;
    double f1 = e1;
    double f2 = e2 + e3;
    double f3 = e3 - e2;
    double f4 = e4;
    double f5 = e5 + e7;
    double f6 = e6;
    double f7 = e7 - e5;
    
    double g0 = f0;
    double g1 = f1;
    double g2 = f2;
    double g3 = f3;
    double g4 = f4 + f7;
    double g5 = f5 + f6;
    double g6 = -f6 + f5;
    double g7 = f7 - f4;
    
#if 1
    const double s0 = (cos(0)*sqrt(.5)/2)/(1       );  // 0.353553
    const double s1 = (cos(1.*M_PI/16)/2)/(-a5+a4+1);  // 0.254898
    const double s2 = (cos(2.*M_PI/16)/2)/(a1+1    );  // 0.270598
    const double s3 = (cos(3.*M_PI/16)/2)/(a5+1    );  // 0.300672
    const double s4 = s0;  // (cos(4.*M_PI/16)/2)/(1       );
    const double s5 = (cos(5.*M_PI/16)/2)/(1-a5    );  // 0.449988
    const double s6 = (cos(6.*M_PI/16)/2)/(1-a1    );  // 0.653281
    const double s7 = (cos(7.*M_PI/16)/2)/(a5-a4+1 );  // 1.281458
#else
    const double s0 = 1.;
    const double s1 = 1.;
    const double s2 = 1.;
    const double s3 = 1.;
    const double s4 = 1.;
    const double s5 = 1.;
    const double s6 = 1.;
    const double s7 = 1.;
#endif
    
    o[0] = g0 * s0;
    o[4] = g1 * s4;
    o[2] = g2 * s2;
    o[6] = g3 * s6;
    o[5] = g4 * s5;
    o[1] = g5 * s1;
    o[7] = g6 * s7;
    o[3] = g7 * s3;
}


void llm_dct_don(const double in[8], double out[8]) {
    //constants:
    const double sin1 = sin(1. * M_PI / 16.);
    const double cos1 = cos(1. * M_PI / 16.);
    const double sin3 = sin(3. * M_PI / 16.);
    const double cos3 = cos(3. * M_PI / 16.);
    const double root2sin6 = sqrt(2.) * sin(6. * M_PI / 16.);
    const double root2cos6 = sqrt(2.) * cos(6. * M_PI / 16.);
    
    //after stage 1:
    const double sum1_0 = in[0] + in[7];
    const double sum1_1 = in[1] + in[6];
    const double sum1_2 = in[2] + in[5];
    const double sum1_3 = in[3] + in[4];
    const double sum1_4 = -in[4] + in[3];
    const double sum1_5 = -in[5] + in[2];
    const double sum1_6 = -in[6] + in[1];
    const double sum1_7 = -in[7] + in[0];
    
    //after stage 2:
    const double sum2_0 = sum1_0 + sum1_3;
    const double sum2_1 = sum1_1 + sum1_2;
    const double sum2_2 = -sum1_2 + sum1_1;
    const double sum2_3 = -sum1_3 + sum1_0;
    
    const double z1 = cos3 * (sum1_7 + sum1_4);
    const double sum2_4 = (sin3 - cos3) * sum1_7 + z1;
    const double sum2_7 = (-sin3 - cos3) * sum1_4 + z1;
    
    const double z2 = cos1 * (sum1_6 + sum1_5);
    const double sum2_5 = (sin1 - cos1) * sum1_6 + z2;
    const double sum2_6 = (-sin1 - cos1) * sum1_5 + z2;
    
    //after stage 3:
    const double sum3_0 = sum2_0 + sum2_1;
    const double sum3_1 = -sum2_1 + sum2_0;
    
    const double z3 = root2cos6 * (sum2_3 + sum2_2);
    const double sum3_2 = (root2cos6 - root2cos6) * sum2_3 + z3;
    const double sum3_3 = (-root2cos6 - root2cos6) * sum2_2 + z3;
    
    const double sum3_4 = sum2_4 + sum2_6;
    const double sum3_5 = -sum2_5 + sum2_7;
    const double sum3_6 = -sum2_6 + sum2_4;
    const double sum3_7 = sum2_7 + sum2_5;
    
    //after stage 4:
    const double sum4_4 = -sum3_4 + sum3_7;
    const double sum4_5 = sum3_5 * sqrt(2.);
    const double sum4_6 = sum3_6 * sqrt(2.);
    const double sum4_7 = sum3_7 + sum3_4;
    
    //shuffle and scaling:
    out[0] = sum3_0 / sqrt(8.);
    out[4] = sum3_1 / sqrt(8.);
    out[2] = sum3_2 / sqrt(8.);
    out[6] = sum3_3 / sqrt(8.);
    out[7] = sum4_4 / sqrt(8.);
    out[3] = sum4_5 / sqrt(8.);
    out[5] = sum4_6 / sqrt(8.);
    out[1] = sum4_7 / sqrt(8.);
    
    
    
    
}