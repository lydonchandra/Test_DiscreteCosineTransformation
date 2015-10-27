//
//  dct_2d.c
//  Test_DiscreteCosineTransformation
//
//  Created by Don on 10/23/15.
//  Copyright © 2015 Don. All rights reserved.
//

#include "dct_2d.h"

/* DCT and IDCT - listing 3
 * Copyright (c) 2001 Emil Mikulic.
 * http://unix4lyfe.org/dct/
 *
 * Feel free to do whatever you like with this code.
 * Feel free to credit me.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "targa.h"



typedef uint8_t byte;

#define pixel(i,x,y) ( (i)->image_data[((y)*( (i)->width ))+(x)] )

#define DONTFAIL(x) do { tga_result res;  if ((res=x) != TGA_NOERR) { \
printf("Targa error: %s\n", tga_error(res)); \
exit(EXIT_FAILURE); } } while(0)



void load_tga(tga_image *tga, const char *fn)
{
    DONTFAIL( tga_read(tga, fn) );
    
//    printf("Loaded %dx%dx%dbpp targa (\"%s\").\n",
//           tga->width, tga->height, tga->pixel_depth, fn);
    
    if (!tga_is_mono(tga)) DONTFAIL( tga_desaturate_rec_601_1(tga) );
    if (!tga_is_top_to_bottom(tga)) DONTFAIL( tga_flip_vert(tga) );
    if (tga_is_right_to_left(tga)) DONTFAIL( tga_flip_horiz(tga) );
    
    if ((tga->width % 8 != 0) || (tga->height % 8 != 0))
    {
//        printf("Width and height must be multiples of 8\n");
        exit(EXIT_FAILURE);
    }
}



#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265358979
#endif
#endif



/* Fast DCT algorithm due to Arai, Agui, Nakajima
 * Implementation due to Tim Kientzle
 */
void dct(tga_image *tga, double data[8][8],
         const int xpos, const int ypos)
{
    int i;
    int rows[8][8];
    
    static const int	c1=1004 /* cos(pi/16) << 10 */,
				s1=200 /* sin(pi/16) */,
				c3=851 /* cos(3pi/16) << 10 */,
				s3=569 /* sin(3pi/16) << 10 */,
				r2c6=554 /* sqrt(2)*cos(6pi/16) << 10 */,
				r2s6=1337 /* sqrt(2)*sin(6pi/16) << 10 */,
				r2=181; /* sqrt(2) << 7*/
    
    int x0,x1,x2,x3,x4,x5,x6,x7,x8;
    
    /* transform rows */
    for (i=0; i<8; i++)
    {
        x0 = pixel(tga, xpos+0, ypos+i);
        x1 = pixel(tga, xpos+1, ypos+i);
        x2 = pixel(tga, xpos+2, ypos+i);
        x3 = pixel(tga, xpos+3, ypos+i);
        x4 = pixel(tga, xpos+4, ypos+i);
        x5 = pixel(tga, xpos+5, ypos+i);
        x6 = pixel(tga, xpos+6, ypos+i);
        x7 = pixel(tga, xpos+7, ypos+i);
        
        /* Stage 1 */
        x8=x7+x0;
        x0-=x7;
        x7=x1+x6;
        x1-=x6;
        x6=x2+x5;
        x2-=x5;
        x5=x3+x4;
        x3-=x4;
        
        /* Stage 2 */
        x4=x8+x5;
        x8-=x5;
        x5=x7+x6;
        x7-=x6;
        x6=c1*(x1+x2);
        x2=(-s1-c1)*x2+x6;
        x1=(s1-c1)*x1+x6;
        x6=c3*(x0+x3);
        x3=(-s3-c3)*x3+x6;
        x0=(s3-c3)*x0+x6;
        
        /* Stage 3 */
        x6=x4+x5;
        x4-=x5;
        x5=r2c6*(x7+x8);
        x7=(-r2s6-r2c6)*x7+x5;
        x8=(r2s6-r2c6)*x8+x5;
        x5=x0+x2;
        x0-=x2;
        x2=x3+x1;
        x3-=x1;
        
        /* Stage 4 and output */
        rows[i][0]=x6;
        rows[i][4]=x4;
        rows[i][2]=x8>>10;
        rows[i][6]=x7>>10;
        rows[i][7]=(x2-x5)>>10;
        rows[i][1]=(x2+x5)>>10;
        rows[i][3]=(x3*r2)>>17;
        rows[i][5]=(x0*r2)>>17;
    }
    
    /* transform columns */
    for (i=0; i<8; i++)
    {
        x0 = rows[0][i];
        x1 = rows[1][i];
        x2 = rows[2][i];
        x3 = rows[3][i];
        x4 = rows[4][i];
        x5 = rows[5][i];
        x6 = rows[6][i];
        x7 = rows[7][i];
        
        /* Stage 1 */
        x8=x7+x0;
        x0-=x7;
        x7=x1+x6;
        x1-=x6;
        x6=x2+x5;
        x2-=x5;
        x5=x3+x4;
        x3-=x4;
        
        /* Stage 2 */
        x4=x8+x5;
        x8-=x5;
        x5=x7+x6;
        x7-=x6;
        x6=c1*(x1+x2);
        x2=(-s1-c1)*x2+x6;
        x1=(s1-c1)*x1+x6;
        x6=c3*(x0+x3);
        x3=(-s3-c3)*x3+x6;
        x0=(s3-c3)*x0+x6;
        
        /* Stage 3 */
        x6=x4+x5;
        x4-=x5;
        x5=r2c6*(x7+x8);
        x7=(-r2s6-r2c6)*x7+x5;
        x8=(r2s6-r2c6)*x8+x5;
        x5=x0+x2;
        x0-=x2;
        x2=x3+x1;
        x3-=x1;
        
        /* Stage 4 and output */
        data[0][i]=(double)((x6+16)>>3);
        data[4][i]=(double)((x4+16)>>3);
        data[2][i]=(double)((x8+16384)>>13);
        data[6][i]=(double)((x7+16384)>>13);
        data[7][i]=(double)((x2-x5+16384)>>13);
        data[1][i]=(double)((x2+x5+16384)>>13);
        data[3][i]=(double)(((x3>>8)*r2+8192)>>12);
        data[5][i]=(double)(((x0>>8)*r2+8192)>>12);
    }
}



/* play with this bit */
void quantize(double dct_buf[8][8])
{
    int x,y;
    
    for (y=0; y<8; y++)
        for (x=0; x<8; x++)
            if (x > 1 || y > 1) dct_buf[y][x] = 0.0;
}



#define COEFFS(Cu,Cv,u,v) { \
if (u == 0) Cu = 1.0 / sqrt(2.0); else Cu = 1.0; \
if (v == 0) Cv = 1.0 / sqrt(2.0); else Cv = 1.0; \
}

void idct(tga_image *tga, double data[8][8], const int xpos, const int ypos)
{
    int u,v,x,y;
    
    /* iDCT */
    for (y=0; y<8; y++)
        for (x=0; x<8; x++)
        {
            double z = 0.0;
            
            for (v=0; v<8; v++)
                for (u=0; u<8; u++)
                {
                    double S, q;
                    double Cu, Cv;
                    
                    COEFFS(Cu,Cv,u,v);
                    S = data[v][u];
                    
                    q = Cu * Cv * S *
                    cos((double)(2*x+1) * (double)u * PI/16.0) *
                    cos((double)(2*y+1) * (double)v * PI/16.0);
                    
                    z += q;
                }
            
            z /= 4.0;
            if (z > 255.0) z = 255.0;
            if (z < 0) z = 0.0;
            
            pixel(tga, x+xpos, y+ypos) = (uint8_t) z;
        }
}

#define PI_DIV_16 0.19634954084936207740391521145497

#define _2x0_Plus1_x_PI_DIV_16_u0  0.0
#define _2x0_Plus1_x_PI_DIV_16_u1  0.19634954084936207740391521145497
#define _2x0_Plus1_x_PI_DIV_16_u2  0.39269908169872415480783042290994
#define _2x0_Plus1_x_PI_DIV_16_u3  0.58904862254809
#define _2x0_Plus1_x_PI_DIV_16_u4  0.78539816339745
#define _2x0_Plus1_x_PI_DIV_16_u5  0.98174770424681
#define _2x0_Plus1_x_PI_DIV_16_u6  1.17809724509617
#define _2x0_Plus1_x_PI_DIV_16_u7  1.37444678594553
#define cos_2x0_Plus1_x_PI_DIV_16_u0   1.0
#define cos_2x0_Plus1_x_PI_DIV_16_u1   0.98078528
#define cos_2x0_Plus1_x_PI_DIV_16_u2   0.923879533
#define cos_2x0_Plus1_x_PI_DIV_16_u3   0.831469612
#define cos_2x0_Plus1_x_PI_DIV_16_u4   0.707106781
#define cos_2x0_Plus1_x_PI_DIV_16_u5   0.555570233
#define cos_2x0_Plus1_x_PI_DIV_16_u6   0.382683432
#define cos_2x0_Plus1_x_PI_DIV_16_u7   0.195090322

#define _2x1_Plus1_x_PI_DIV_16_u1  0.58904862254809
#define _2x1_Plus1_x_PI_DIV_16_u2  1.17809724509618
#define _2x1_Plus1_x_PI_DIV_16_u3  1.76714586764427
#define _2x1_Plus1_x_PI_DIV_16_u4  2.35619449019236
#define _2x1_Plus1_x_PI_DIV_16_u5  2.94524311274045
#define _2x1_Plus1_x_PI_DIV_16_u6  3.53429173528854
#define _2x1_Plus1_x_PI_DIV_16_u7  4.12334035783663
#define cos_2x1_Plus1_x_PI_DIV_16_u1   0.831469612
#define cos_2x1_Plus1_x_PI_DIV_16_u2   0.382683432
#define cos_2x1_Plus1_x_PI_DIV_16_u3   -0.195090322
#define cos_2x1_Plus1_x_PI_DIV_16_u4   -0.707106781
#define cos_2x1_Plus1_x_PI_DIV_16_u5   -0.98078528
#define cos_2x1_Plus1_x_PI_DIV_16_u6   -0.923879533
#define cos_2x1_Plus1_x_PI_DIV_16_u7   -0.555570233


#define _2x2_Plus1_x_PI_DIV_16_u1  0.98174770424681
#define _2x2_Plus1_x_PI_DIV_16_u2  1.96349540849362
#define _2x2_Plus1_x_PI_DIV_16_u3  2.94524311274043
#define _2x2_Plus1_x_PI_DIV_16_u4  3.92699081698724
#define _2x2_Plus1_x_PI_DIV_16_u5  4.90873852123405
#define _2x2_Plus1_x_PI_DIV_16_u6  5.89048622548086
#define _2x2_Plus1_x_PI_DIV_16_u7  6.87223392972767
#define cos_2x2_Plus1_x_PI_DIV_16_u1 0.555570233
#define cos_2x2_Plus1_x_PI_DIV_16_u2 -0.382683432
#define cos_2x2_Plus1_x_PI_DIV_16_u3 -0.98078528
#define cos_2x2_Plus1_x_PI_DIV_16_u4 -0.707106781
#define cos_2x2_Plus1_x_PI_DIV_16_u5 0.195090322
#define cos_2x2_Plus1_x_PI_DIV_16_u6 0.923879533
#define cos_2x2_Plus1_x_PI_DIV_16_u7 0.831469612


#define _2x3_Plus1_x_PI_DIV_16_u1  1.374446786
#define _2x3_Plus1_x_PI_DIV_16_u2  2.748893572
#define _2x3_Plus1_x_PI_DIV_16_u3  4.123340358
#define _2x3_Plus1_x_PI_DIV_16_u4  5.497787144
#define _2x3_Plus1_x_PI_DIV_16_u5  6.87223393
#define _2x3_Plus1_x_PI_DIV_16_u6  8.246680716
#define _2x3_Plus1_x_PI_DIV_16_u7  9.621127502
#define cos_2x3_Plus1_x_PI_DIV_16_u1    0.195090322
#define cos_2x3_Plus1_x_PI_DIV_16_u2    -0.923879533
#define cos_2x3_Plus1_x_PI_DIV_16_u3    -0.555570233
#define cos_2x3_Plus1_x_PI_DIV_16_u4    0.707106781
#define cos_2x3_Plus1_x_PI_DIV_16_u5    0.831469612
#define cos_2x3_Plus1_x_PI_DIV_16_u6    -0.382683433
#define cos_2x3_Plus1_x_PI_DIV_16_u7    -0.98078528


#define _2x4_Plus1_x_PI_DIV_16_u1  1.76714586764426
#define _2x4_Plus1_x_PI_DIV_16_u2  3.53429173528852
#define _2x4_Plus1_x_PI_DIV_16_u3  5.30143760293278
#define _2x4_Plus1_x_PI_DIV_16_u4  7.06858347057704
#define _2x4_Plus1_x_PI_DIV_16_u5  8.8357293382213
#define _2x4_Plus1_x_PI_DIV_16_u6  10.60287520586556
#define _2x4_Plus1_x_PI_DIV_16_u7  12.37002107350982
#define cos_2x4_Plus1_x_PI_DIV_16_u1    -0.195090322
#define cos_2x4_Plus1_x_PI_DIV_16_u2    -0.923879533
#define cos_2x4_Plus1_x_PI_DIV_16_u3    0.555570233
#define cos_2x4_Plus1_x_PI_DIV_16_u4    0.707106781
#define cos_2x4_Plus1_x_PI_DIV_16_u5    -0.831469612
#define cos_2x4_Plus1_x_PI_DIV_16_u6    -0.382683432
#define cos_2x4_Plus1_x_PI_DIV_16_u7    0.98078528


#define _2x5_Plus1_x_PI_DIV_16_u1  2.15984494934298
#define _2x5_Plus1_x_PI_DIV_16_u2  4.31968989868596
#define _2x5_Plus1_x_PI_DIV_16_u3  6.47953484802894
#define _2x5_Plus1_x_PI_DIV_16_u4  8.63937979737192
#define _2x5_Plus1_x_PI_DIV_16_u5  10.7992247467149
#define _2x5_Plus1_x_PI_DIV_16_u6  12.95906969605788
#define _2x5_Plus1_x_PI_DIV_16_u7  15.11891464540086
#define cos_2x5_Plus1_x_PI_DIV_16_u1   -0.555570233
#define cos_2x5_Plus1_x_PI_DIV_16_u2   -0.382683432
#define cos_2x5_Plus1_x_PI_DIV_16_u3   0.98078528
#define cos_2x5_Plus1_x_PI_DIV_16_u4   -0.707106781
#define cos_2x5_Plus1_x_PI_DIV_16_u5   -0.195090322
#define cos_2x5_Plus1_x_PI_DIV_16_u6   0.923879533
#define cos_2x5_Plus1_x_PI_DIV_16_u7    -0.831469612


#define _2x6_Plus1_x_PI_DIV_16_u1  2.552544031
#define _2x6_Plus1_x_PI_DIV_16_u2  5.105088062
#define _2x6_Plus1_x_PI_DIV_16_u3  7.657632093
#define _2x6_Plus1_x_PI_DIV_16_u4  10.21017612
#define _2x6_Plus1_x_PI_DIV_16_u5  12.76272016
#define _2x6_Plus1_x_PI_DIV_16_u6  15.31526419
#define _2x6_Plus1_x_PI_DIV_16_u7  17.86780822
#define cos_2x6_Plus1_x_PI_DIV_16_u1    -0.831469612
#define cos_2x6_Plus1_x_PI_DIV_16_u2    0.382683432
#define cos_2x6_Plus1_x_PI_DIV_16_u3    0.195090322
#define cos_2x6_Plus1_x_PI_DIV_16_u4    -0.707106784
#define cos_2x6_Plus1_x_PI_DIV_16_u5    0.980785279
#define cos_2x6_Plus1_x_PI_DIV_16_u6    -0.923879534
#define cos_2x6_Plus1_x_PI_DIV_16_u7    0.555570235


#define _2x7_Plus1_x_PI_DIV_16_u1  2.945243113
#define _2x7_Plus1_x_PI_DIV_16_u2  5.890486225
#define _2x7_Plus1_x_PI_DIV_16_u3  8.835729338
#define _2x7_Plus1_x_PI_DIV_16_u4  11.78097245
#define _2x7_Plus1_x_PI_DIV_16_u5  14.72621556
#define _2x7_Plus1_x_PI_DIV_16_u6  17.67145868
#define _2x7_Plus1_x_PI_DIV_16_u7  20.61670179
#define cos_2x7_Plus1_x_PI_DIV_16_u1    -0.98078528
#define cos_2x7_Plus1_x_PI_DIV_16_u2    0.923879532
#define cos_2x7_Plus1_x_PI_DIV_16_u3    -0.831469612
#define cos_2x7_Plus1_x_PI_DIV_16_u4    0.707106781
#define cos_2x7_Plus1_x_PI_DIV_16_u5    -0.55557023
#define cos_2x7_Plus1_x_PI_DIV_16_u6    0.382683436
#define cos_2x7_Plus1_x_PI_DIV_16_u7   -0.195090323




void idct_don(tga_image *tga, double data[8][8], const int xpos, const int ypos)
{
    double cosData[8][8] = {
        {1.0, 0.98078528, 0.923879533, 0.831469612, 0.707106781, 0.555570233, 0.382683432,0.195090322},
        
        {1.0, 0.831469612,0.382683432,-0.195090322,-0.707106781,-0.98078528,-0.923879533,-0.555570233},
        
        {1.0, 0.555570233,-0.382683432,-0.98078528,-0.707106781,0.195090322,0.923879533,0.831469612},
        
        {1.0, 0.195090322,-0.923879533,-0.555570233,0.707106781,0.831469612,-0.382683433,-0.98078528},
        
        {1.0, -0.195090322,-0.923879533,0.555570233,0.707106781,-0.831469612,-0.382683432,0.98078528},
        
        {1.0, -0.555570233,-0.382683432,0.98078528,-0.707106781,-0.195090322,0.923879533,-0.831469612},
        
        {1.0, -0.831469612,0.382683432,0.195090322,-0.707106784,0.980785279,-0.923879534,0.555570235},
        
        {1.0, -0.98078528,0.923879532,-0.831469612,0.707106781,-0.55557023,0.382683436,-0.195090323}
    };
    
    int /*u,*/v,x,y;
    
    uint8_t *image_data = tga->image_data;
    uint16_t image_width = tga->width;
    
    /* iDCT */
    for (y=0; y<8; y++)
    {
        //why using left shift gives us blurry pic?
        //double _2yPlus1_x_PI_DIV_16 = ((y << 2) + 1) * PI_DIV_16;
        double *cosDataY = cosData[y];
        
        int yIndex = ((y+ypos) * image_width);
        uint8_t *image_data_y = &image_data[yIndex];
        
        for (x=0; x<8; x++)
        {
            double z = 0.0;
            double *cosDataX = cosData[x];
            double cosDataX_1 = cosDataX[1];
            double cosDataX_2 = cosDataX[2];
            double cosDataX_3 = cosDataX[3];
            double cosDataX_4 = cosDataX[4];
            double cosDataX_5 = cosDataX[5];
            double cosDataX_6 = cosDataX[6];
            double cosDataX_7 = cosDataX[7];
            
//            for (v=0; v<8; v++)
//            {
            
            //v = 0
            double *data_v = data[0];
            z = z
            +
            (0.70710678118655 * cosDataY[0])
            *
            (   //u = 0
                (data_v[0] * 0.70710678118655) +
                //u = 1;
                (data_v[1] * cosDataX_1) +
                //u = 2;
                (data_v[2] * cosDataX_2) +
                //u = 3;
                (data_v[3] * cosDataX_3) +
                //u = 4;
                (data_v[4] * cosDataX_4) +
                //u = 5;
                (data_v[5] * cosDataX_5) +
                //u = 6;
                (data_v[6] * cosDataX_6) +
                //u = 7;
                (data_v[7] * cosDataX_7)
             );
            
            //v = 1
            data_v = data[1];
            z = z
            +
            cosDataY[1]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );
            
            //v = 2
            data_v = data[2];
            z = z
            +
            cosDataY[2]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );
            
            //v = 3
            data_v = data[3];
            z = z
            +
            cosDataY[3]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );
            
            //v = 4
            data_v = data[4];
            z = z
            +
            cosDataY[4]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );
            
            //v = 5
            data_v = data[5];
            z = z
            +
            cosDataY[5]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );
            
            //v = 6
            data_v = data[6];
            z = z
            +
            cosDataY[6]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );
            
            //v = 7
            data_v = data[7];
            z = z
            +
            cosDataY[7]
            *
            (   //u = 0
             (data_v[0] * 0.70710678118655) +
             //u = 1;
             (data_v[1] * cosDataX_1) +
             //u = 2;
             (data_v[2] * cosDataX_2) +
             //u = 3;
             (data_v[3] * cosDataX_3) +
             //u = 4;
             (data_v[4] * cosDataX_4) +
             //u = 5;
             (data_v[5] * cosDataX_5) +
             //u = 6;
             (data_v[6] * cosDataX_6) +
             //u = 7;
             (data_v[7] * cosDataX_7)
             );

//            double /*Cu = 1.0,*/ Cv = 1.0;
//            double cos_2yPlus1_x_PI_DIV_16_x_v = cosDataY[v];
//            double *data_v = data[v];
//            
//            if (v == 0) Cv = 0.70710678118655; else Cv = 1.0;
//            double Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v = Cv * cos_2yPlus1_x_PI_DIV_16_x_v;
//            
//            z = z
//            +
//            Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v
//            *
//            (   //u = 0
//             (data_v[0] * 0.70710678118655)
//             +
//             //u = 1;
//             (data_v[1] * cosDataX[1])
//             +
//             //u = 2;
//             (data_v[2] * cosDataX[2])
//             +
//             //u = 3;
//             (data_v[3] * cosDataX[3])
//             +
//             //u = 4;
//             (data_v[4] * cosDataX[4])
//             +
//             //u = 5;
//             (data_v[5] * cosDataX[5])
//             +
//             //u = 6;
//             (data_v[6] * cosDataX[6])
//             +
//             //u = 7;
//             (data_v[7] *cosDataX[7])
//             );

            
            
            
//            } //end for(v...
            
            z /= 4.0;
            if (z > 255.0) z = 255.0;
            if (z < 0) z = 0.0;
            
            //pixel(tga, x+xpos, y+ypos) = (uint8_t) z;
            //image_data[ yIndex + (x+xpos) ] = (uint8_t)z;
            image_data_y[ (x+xpos) ] = (uint8_t)z;
            
            //            #define pixel(i,x,y) ( (i)->image_data[((y)*( (i)->width ))+(x)] )
        }
    }
}




void idct_don2(tga_image *tga, double data[8][8], const int xpos, const int ypos)
{
    double cosData[8][8] = {
        {1.0, 0.98078528, 0.923879533, 0.831469612, 0.707106781, 0.555570233, 0.382683432,0.195090322},
        
        {1.0, 0.831469612,0.382683432,-0.195090322,-0.707106781,-0.98078528,-0.923879533,-0.555570233},
        
        {1.0, 0.555570233,-0.382683432,-0.98078528,-0.707106781,0.195090322,0.923879533,0.831469612},
        
        {1.0, 0.195090322,-0.923879533,-0.555570233,0.707106781,0.831469612,-0.382683433,-0.98078528},
        
        {1.0, -0.195090322,-0.923879533,0.555570233,0.707106781,-0.831469612,-0.382683432,0.98078528},
        
        {1.0, -0.555570233,-0.382683432,0.98078528,-0.707106781,-0.195090322,0.923879533,-0.831469612},
        
        {1.0, -0.831469612,0.382683432,0.195090322,-0.707106784,0.980785279,-0.923879534,0.555570235},
        
        {1.0, -0.98078528,0.923879532,-0.831469612,0.707106781,-0.55557023,0.382683436,-0.195090323}
    };
    
    int /*u,*/v,x,y;
    
    /* iDCT */
    for (y=0; y<8; y++)
    {
        //double _2yPlus1_x_PI_DIV_16 = (2. * y + 1) * PI_DIV_16;
        //why using left shift gives us blurry pic?
        //double _2yPlus1_x_PI_DIV_16 = ((y << 2) + 1) * PI_DIV_16;
        double *cosDataY = cosData[y];
        
        for (x=0; x<8; x++)
        {
            double z = 0.0;
            
//            double _2xPlus1_x_PI_DIV_16 = (2. * x + 1) * PI_DIV_16;
//            double cos_2xPlus1_x_PI_DIV_16_x1 = cos(_2xPlus1_x_PI_DIV_16);
//            double cos_2xPlus1_x_PI_DIV_16_x2 = cos(_2xPlus1_x_PI_DIV_16 * 2.);
//            double cos_2xPlus1_x_PI_DIV_16_x3 = cos(_2xPlus1_x_PI_DIV_16 * 3.);
//            double cos_2xPlus1_x_PI_DIV_16_x4 = cos(_2xPlus1_x_PI_DIV_16 * 4.);
//            double cos_2xPlus1_x_PI_DIV_16_x5 = cos(_2xPlus1_x_PI_DIV_16 * 5.);
//            double cos_2xPlus1_x_PI_DIV_16_x6 = cos(_2xPlus1_x_PI_DIV_16 * 6.);
//            double cos_2xPlus1_x_PI_DIV_16_x7 = cos(_2xPlus1_x_PI_DIV_16 * 7.);
            double *cosDataX = cosData[x];
            

            for (v=0; v<8; v++)
                //for (u=0; u<8; u++)
                {
//#define COEFFS(Cu,Cv,u,v) { \
//if (u == 0) Cu = 1.0 / sqrt(2.0); else Cu = 1.0; \
//if (v == 0) Cv = 1.0 / sqrt(2.0); else Cv = 1.0; \
//}
                    double S, q;
                    double /*Cu = 1.0,*/ Cv = 1.0;
                    
                    //double cos_2yPlus1_x_PI_DIV_16_x_v = cos(_2yPlus1_x_PI_DIV_16 * (double)v);
                    double cos_2yPlus1_x_PI_DIV_16_x_v = cosDataY[v];
                    double *data_v = data[v];
                    
//                    u = 0;
                    //COEFFS(Cu,Cv,u,v);
                    //Cu = 0.70710678118655;
                    if (v == 0) Cv = 0.70710678118655; else Cv = 1.0;
                    double Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v = Cv * cos_2yPlus1_x_PI_DIV_16_x_v;
                    //Cv = 0.70710678118655;
                    S = data_v[0];
                    
                    q = 0.70710678118655 * Cv * S *
                    //1.
                    //*
                    cos_2yPlus1_x_PI_DIV_16_x_v;
                    
                    z += q;
                    
//                    q = Cu * Cv * S *
//                    cos((double)(2*x+1) * (double)u * PI/16.0) *
//                    cos((double)(2*y+1) * (double)v * PI/16.0);
                    
                    //Cu = Cv = 1.0;
                    //u = 1;
                    //COEFFS(Cu, Cv, 1, v);
                    //Cu = 1.0;
                    //double Cu_x_Cv = /*1 **/ Cv;
                    //S = data_v[1];
                    z += data_v[1] *
                        cosDataX[1]
//                        cos_2xPlus1_x_PI_DIV_16_x1
                        *
                        Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;
                    
                    //u = 2;
                    //COEFFS(Cu, Cv, u, v);
                    //S = data_v[2];
                    z += data_v[2] *
                        cosDataX[2]
//                    cos_2xPlus1_x_PI_DIV_16_x2
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;

                    //u = 3;
                    //COEFFS(Cu, Cv, u, v);
                    //S = data_v[3];
                    z += data_v[3] *
                        cosDataX[3]
//                    cos_2xPlus1_x_PI_DIV_16_x3
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;

                    //u = 4;
                    //COEFFS(Cu, Cv, u, v);
                    //S = data_v[4];
                    z += data_v[4] *
                        cosDataX[4]
//                    cos_2xPlus1_x_PI_DIV_16_x4
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;

                    //u = 5;
                    //COEFFS(Cu, Cv, u, v);
                    //S = data_v[5];
                    z += data_v[5] *
                    cosDataX[5]
//                    cos_2xPlus1_x_PI_DIV_16_x5
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;

                    //u = 6;
                    //COEFFS(Cu, Cv, u, v);
                    //S = data_v[6];
                    z += data_v[6] *
                    cosDataX[6]
//                    cos_2xPlus1_x_PI_DIV_16_x6
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;

                    //u = 7;
                    //COEFFS(Cu, Cv, u, v);
                    //S = data_v[7];
                    z += data_v[7] *
                    cosDataX[7]
//                    cos_2xPlus1_x_PI_DIV_16_x7
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    //z += q;

                } //end for(v...
            
            z /= 4.0;
            if (z > 255.0) z = 255.0;
            if (z < 0) z = 0.0;
            
            //pixel(tga, x+xpos, y+ypos) = (uint8_t) z;
            (tga)->image_data[ ((y+ypos) * (tga)->width) + (x+xpos) ] = (uint8_t)z;
            
//            #define pixel(i,x,y) ( (i)->image_data[((y)*( (i)->width ))+(x)] )
        }
    }
}



int convert_to_tga_mono(const char* filePath, const char* outFilePath)
{
    tga_image tga;
    double dct_buf[8][8];
    int i, j, k, l;
    
    load_tga(&tga, filePath);
    k = 0;
    l = (tga.height / 8) * (tga.width / 8);
    for (j=0; j<tga.height/8; j++) {
        for (i=0; i<tga.width; i++) {
            dct(&tga, dct_buf, i*8, j*8);
            quantize(dct_buf);
            idct(&tga, dct_buf, i*8, j*8);
        }
    }
    
    tga_write_mono(outFilePath, tga.image_data, tga.width, tga.height);
    tga_free_buffers(&tga);
    return EXIT_SUCCESS;
}


int main_old2(const char* filePath, const char* outFilePath)
{
    tga_image tga;
    double dct_buf[8][8];
    int i, j, k, l;
    
//    load_tga(&tga, "in.tga");
    load_tga(&tga, filePath);
    k = 0;
    l = (tga.height / 8) * (tga.width / 8);
    for (j=0; j<tga.height/8; j++)
        for (i=0; i<tga.width/8; i++)
        {
            dct(&tga, dct_buf, i*8, j*8);
            quantize(dct_buf);
            //idct(&tga, dct_buf, i*8, j*8);
            idct_don(&tga, dct_buf, i*8, j*8);
//            printf("processed %d/%d blocks.\r", ++k,l);
//            fflush(stdout);
        }
//    printf("\n");
    
    DONTFAIL( tga_write_mono(outFilePath, tga.image_data,
                             tga.width, tga.height) );
    
    tga_free_buffers(&tga);
    return EXIT_SUCCESS;
}


//idct -Ofast disassembly code
//0x1000ac6d0 <+0>:   movz   x8, #0
//0x1000ac6d4 <+4>:   fmov   d0, #1.00000000
//0x1000ac6d8 <+8>:   nop
//0x1000ac6dc <+12>:  ldr    d1, 0x1000ae7c0
//0x1000ac6e0 <+16>:  adr    x9, #8424
//0x1000ac6e4 <+20>:  nop
//0x1000ac6e8 <+24>:  fmov   d2, #0.25000000
//->  0x1000ac6ec <+28>:  ubfx   x10, x2, #0, #32
//0x1000ac6f0 <+32>:  ubfx   x11, x3, #0, #32
//0x1000ac6f4 <+36>:  add    x12, x1, #32
//0x1000ac6f8 <+40>:  nop
//0x1000ac6fc <+44>:  ldr    d3, 0x1000ae7b8
//0x1000ac700 <+48>:  mov    x13, x9
//0x1000ac704 <+52>:  movz   x14, #0
//0x1000ac708 <+56>:  add    x15, x8, x11
//0x1000ac70c <+60>:  movz   x16, #0
//0x1000ac710 <+64>:  add    x1, x9, x14, lsl #6
//0x1000ac714 <+68>:  ldp    d4, d5, [x1, #8]
//0x1000ac718 <+72>:  ldp    d6, d7, [x1, #24]
//0x1000ac71c <+76>:  ldp    d16, d17, [x1, #40]
//0x1000ac720 <+80>:  movi.16b v18, #0
//0x1000ac724 <+84>:  mov    x17, x13
//0x1000ac728 <+88>:  ldr    d19, [x1, #56]
//0x1000ac72c <+92>:  ldr    d20, [x17], #8
//0x1000ac730 <+96>:  cmp    w16, #0
//0x1000ac734 <+100>: fcsel  d21, d1, d0, eq
//0x1000ac738 <+104>: fmul   d22, d20, d21
//0x1000ac73c <+108>: add    x1, x12, x16
//0x1000ac740 <+112>: fmul   d21, d21, d1
//0x1000ac744 <+116>: fmul   d20, d21, d20
//0x1000ac748 <+120>: ldp    d23, d21, [x1, #-32]
//0x1000ac74c <+124>: fmul   d21, d4, d21
//0x1000ac750 <+128>: ldp    d24, d25, [x1, #-16]
//0x1000ac754 <+132>: ldp    d26, d27, [x1]
//0x1000ac758 <+136>: ldp    d28, d29, [x1, #16]
//0x1000ac75c <+140>: fmadd  d21, d5, d24, d21
//0x1000ac760 <+144>: fmadd  d21, d6, d25, d21
//0x1000ac764 <+148>: fmadd  d21, d7, d26, d21
//0x1000ac768 <+152>: fmadd  d21, d16, d27, d21
//0x1000ac76c <+156>: fmadd  d21, d17, d28, d21
//0x1000ac770 <+160>: fmadd  d21, d19, d29, d21
//0x1000ac774 <+164>: fmadd  d18, d20, d23, d18
//0x1000ac778 <+168>: fmadd  d18, d22, d21, d18
//0x1000ac77c <+172>: add    x16, x16, #64
//0x1000ac780 <+176>: cmp    x16, #512
//0x1000ac784 <+180>: b.ne   0x1000ac72c               ; <+92> at dct_2d.c:439
//0x1000ac788 <+184>: fmul   d4, d18, d2
//0x1000ac78c <+188>: fmin   d4, d4, d3
//0x1000ac790 <+192>: ldrh   w16, [x0, #14]
//0x1000ac794 <+196>: fcvtzs w17, d4
//0x1000ac798 <+200>: add    w1, w14, w10
//0x1000ac79c <+204>: fcmp   d4, #0.0
//0x1000ac7a0 <+208>: csel   w17, wzr, w17, lt
//0x1000ac7a4 <+212>: madd   w16, w16, w15, w1
//0x1000ac7a8 <+216>: ldr    x1, [x0, #40]
//0x1000ac7ac <+220>: strb   w17, [x1, w16, sxtw]
//0x1000ac7b0 <+224>: add    x14, x14, #1
//0x1000ac7b4 <+228>: cmp    x14, #8
//0x1000ac7b8 <+232>: b.ne   0x1000ac70c               ; <+60> at dct_2d.c:540
//0x1000ac7bc <+236>: add    x8, x8, #1
//0x1000ac7c0 <+240>: add    x13, x13, #64
//0x1000ac7c4 <+244>: cmp    x8, #8
//0x1000ac7c8 <+248>: b.ne   0x1000ac704               ; <+52> at dct_2d.c:406
//0x1000ac7cc <+252>: ret

