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
#define _2x1_Plus1_x_PI_DIV_16_x0  1.19634954084936207740391521145497
#define cos_2x1_Plus1_x_PI_DIV_16_x0   0.36575770316999429021426672670326

#define _2x1_Plus1_x_PI_DIV_16_x1  0.58904862254809
#define cos_2x1_Plus1_x_PI_DIV_16_x1   0.8314696123025452370787883776179

#define _2x2_Plus1_x_PI_DIV_16  0.98174770424681
#define cos_2x2_Plus1_x_PI_DIV_16_x2    -0.998356439317077

void idct_don(tga_image *tga, double data[8][8], const int xpos, const int ypos)
{
    int /*u,*/v,x,y;
    
    /* iDCT */
    for (y=0; y<8; y++)
    {
        double _2yPlus1_x_PI_DIV_16 = (2. * y + 1) * PI_DIV_16;
    
        for (x=0; x<8; x++)
        {
            double z = 0.0;
            
            double _2xPlus1_x_PI_DIV_16 = (2. * x + 1) * PI_DIV_16;
            double cos_2xPlus1_x_PI_DIV_16_x1 = cos(_2xPlus1_x_PI_DIV_16);
            double cos_2xPlus1_x_PI_DIV_16_x2 = cos(_2xPlus1_x_PI_DIV_16 * 2.);
            double cos_2xPlus1_x_PI_DIV_16_x3 = cos(_2xPlus1_x_PI_DIV_16 * 3.);
            double cos_2xPlus1_x_PI_DIV_16_x4 = cos(_2xPlus1_x_PI_DIV_16 * 4.);
            double cos_2xPlus1_x_PI_DIV_16_x5 = cos(_2xPlus1_x_PI_DIV_16 * 5.);
            double cos_2xPlus1_x_PI_DIV_16_x6 = cos(_2xPlus1_x_PI_DIV_16 * 6.);
            double cos_2xPlus1_x_PI_DIV_16_x7 = cos(_2xPlus1_x_PI_DIV_16 * 7.);

            for (v=0; v<8; v++)
                //for (u=0; u<8; u++)
                {
//#define COEFFS(Cu,Cv,u,v) { \
//if (u == 0) Cu = 1.0 / sqrt(2.0); else Cu = 1.0; \
//if (v == 0) Cv = 1.0 / sqrt(2.0); else Cv = 1.0; \
//}
                    double S, q;
                    double /*Cu = 1.0,*/ Cv = 1.0;
                    
                    double cos_2yPlus1_x_PI_DIV_16_x_v = cos(_2yPlus1_x_PI_DIV_16 * (double)v);
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
                    
                    
                    //Cu = Cv = 1.0;
                    //u = 1;
                    //COEFFS(Cu, Cv, 1, v);
                    //Cu = 1.0;
                    //double Cu_x_Cv = /*1 **/ Cv;
                    S = data_v[1];
                    q = S *
                        cos_2xPlus1_x_PI_DIV_16_x1
                        *
                        Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;
                    
                    //u = 2;
                    //COEFFS(Cu, Cv, u, v);
                    S = data_v[2];
                    q = S *
                    cos_2xPlus1_x_PI_DIV_16_x2
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;

                    //u = 3;
                    //COEFFS(Cu, Cv, u, v);
                    S = data_v[3];
                    q = S *
                    cos_2xPlus1_x_PI_DIV_16_x3
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;

                    //u = 4;
                    //COEFFS(Cu, Cv, u, v);
                    S = data_v[4];
                    q = S *
                    cos_2xPlus1_x_PI_DIV_16_x4
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;

                    //u = 5;
                    //COEFFS(Cu, Cv, u, v);
                    S = data_v[5];
                    q = S *
                    cos_2xPlus1_x_PI_DIV_16_x5
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;

                    //u = 6;
                    //COEFFS(Cu, Cv, u, v);
                    S = data_v[6];
                    q = S *
                    cos_2xPlus1_x_PI_DIV_16_x6
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;

                    //u = 7;
                    //COEFFS(Cu, Cv, u, v);
                    S = data_v[7];
                    q = S *
                    cos_2xPlus1_x_PI_DIV_16_x7
                    *
                    Cv_x_cos_2yPlus1_x_PI_DIV_16_x_v;
                    z += q;

                }
            
            z /= 4.0;
            if (z > 255.0) z = 255.0;
            if (z < 0) z = 0.0;
            
            pixel(tga, x+xpos, y+ypos) = (uint8_t) z;
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
