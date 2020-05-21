#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tga.h"
#include "model.h"

void swap(int *a, int *b);
void swapf(double *a, double *b);
int abs(int a);
int sign(int a);

// Triangle rasterization algorithm using a scan line and z-buffer 
void triangle(tgaImage *image,
              int x0, int y0, double z0, double u0, double v0,
              int x1, int y1, double z1, double u1, double v1,
              int x2, int y2, double z2, double u2, double v2, double zbuffer[image->width][image->height],
              double intensity, Model *model);

// Drawing with simple lighting model
void render(tgaImage *image, Model *model);

int main(int argc, char **argv)
{
    int rv = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <objfile> <outfile>\n", argv[0]);
        return -1;
    }

    Model *model = NULL;
    tgaImage *image = NULL;

    do {
        model = loadFromObj(argv[1]);
        if (!model) {
            perror("loadFromObj");
            rv = -1;
            break;
        }

        image = tgaNewImage(1000, 1000, RGB);
        if (!image) {
            perror("tgaNewImage");
            rv = -1;
            break;
        }
        loadDiffuseMap(model, argv[3]);
        // Drawing
        render(image, model);

        if (-1 == tgaSaveToFile(image, argv[2])) {
            perror("tgaSateToFile");
            rv = -1;
            break;
        }
    } while (0);

    if (model) {
        freeModel(model);
    }
    if (image) {
        tgaFreeImage(image);
    } 
    return rv;
}

void triangle(tgaImage *image,
              int x0, int y0, double z0, double u0, double v0,
              int x1, int y1, double z1, double u1, double v1,
              int x2, int y2, double z2, double u2, double v2, double zbuffer[image->width][image->height],
              double intensity, Model *model)
{
    tgaColor color;
    // Sorting vertices by 'y' coord
    if (y0 > y1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
        swapf(&z0, &z1);
        swapf(&u0, &u1);
        swapf(&v0, &v1);
    }
    if (y0 > y2) {
        swap(&x0, &x2);
        swap(&y0, &y2);
        swapf(&z0, &z2);
        swapf(&u0, &u2);
        swapf(&v0, &v2);
    }
    if (y1 > y2) {
        swap(&x1, &x2);
        swap(&y1, &y2);
        swapf(&z1, &z2);
        swapf(&u1, &u2);
        swapf(&v1, &v2);
    }

    double tless, tmore, tz, z, za, zb;
    double xleft, xright, x;
    // For texture
    double ua, ub, va, vb, u, v;
    unsigned int h = model->diffuse_map->height;
    unsigned int w = model->diffuse_map->width;
    int r, g, b;
     
    for(int y = y0; y <= y2; y++){
        // Finding boundary coords using the parametric equation of the line
        tmore = (double)(y - y0) / (y2 - y0);
        if(y > y1){
            tless = (double)(y - y1) / (y2 - y1);
            xleft = x2 * tless + x1 * (1 - tless);
            za = z2 * tless + z1 * (1 - tless);
            ua = u2 * tless + u1 * (1 - tless);
            va = v2 * tless + v1 * (1 - tless);
        }
        else{
            tless = (double)(y - y0) / (y1 - y0);
            xleft = x1 * tless + x0 * (1 - tless);
            za = z1 * tless + z0 * (1 - tless);
            ua = u1 * tless + u0 * (1 - tless);
            va = v1 * tless + v0 * (1 - tless);
        }
        xright = x2 * tmore + x0 * (1 - tmore); 
        zb = z2 * tmore + z0 * (1 - tmore);
        ub = u2 * tmore + u0 * (1 - tmore);
        vb = v2 * tmore + v0 * (1 - tmore);
        if (xleft > xright){
            swapf(&xleft, &xright);
        }
        if (za > zb){
            swapf(&za, &zb);
        }
        if (ua > ub){
            swapf(&ua, &ub);
        }
        if (va > vb){
            swapf(&va, &vb);
        }
        x = xleft;
        // Scan line for x
        while(x <= xright){
            tz = (double)(x - xleft) / (xright - xleft);
            z = zb * tz + za * (1 - tz);
            // Checking z-buffer
            if(z > zbuffer[(int)x][y]){
                u = ub * tz + ua * (1 - tz);
                v = vb * tz + va * (1 - tz);
                zbuffer[(int)x][y] = z;
                color = tgaGetPixel(model->diffuse_map, w * u, h * v);
                r = Red(color) * intensity;
                g = Green(color) * intensity;
                b = Blue(color) * intensity;
                tgaSetPixel(image, x, y, tgaRGB(r, g, b));
            }
            x += 1;
        }
    }
}

void render(tgaImage *image, Model *model)
{
    int face, vert;
    int h = image->height;
    int w = image->width;
    Vec3 *p[3];
    int screen_coords[3][2];
    double a[3], b[3], n[3], norm, intensity;
    //int light[3] = {0, 1, 0};
    int light[3] = {0, 0, -1}; //cat
    Vec3 *s[3];

    // Z-buffer creation
    double zbuffer[w][h];
    for(int i = 0; i < w; i++){
        for (int j = 0; j < h; j++){
            zbuffer[i][j] = -1000000;
        }
    }

    for (face = 0; face < model->nface; ++face) {

        // Translating into screen coordinates
        for (vert = 0; vert < 3; ++vert) {
            p[vert] = getVertex(model, face, vert);
            //screen_coords[vert][0] = ((*p[vert])[0] / 4000  + 1.0) * w / 2;
            //screen_coords[vert][1] = ((*p[vert])[2] / 4000 + 1.0) * h / 2;
            screen_coords[vert][0] = ((*p[vert])[0] + 1.0) * w / 2; //cat
            screen_coords[vert][1] = ((*p[vert])[1] + 1.0) * h / 2; //cat
            s[vert] = getDiffuseUV(model, face, vert);
        }

        // Calculating coords of vectors (a, b) in triangle
        for(int i = 0; i < 3; i++){
            a[i] = (*p[1])[i] - (*p[0])[i];
            b[i] = (*p[2])[i] - (*p[0])[i];
        }

        // Calculating normal coords
        n[0] = a[1] * b[2] - a[2] * b[1];
        n[1] = -(a[0] * b[2] - a[2] * b[0]);
        n[2] = a[0] * b[1] - a[1] * b[0];

        // Computing norm of normal
        norm = sqrt(n[0] * n [0] + n[1] * n [1] + n[2] * n [2]);
        
        // Computing light intensity
        intensity = (light[0] * n[0] + light[1] * n[1] + light[2] * n[2]) / norm;
        if(intensity <= 0){
            // Fill the triangle
            triangle(image, screen_coords[0][0], screen_coords[0][1], (*p[0])[2] / 4000, (*s[0])[0], (*s[0])[1],
                            screen_coords[1][0], screen_coords[1][1], (*p[1])[2] / 4000, (*s[1])[0], (*s[1])[1],
                            screen_coords[2][0], screen_coords[2][1], (*p[2])[2] / 4000, (*s[2])[0], (*s[2])[1], 
                            zbuffer, fabs(intensity), model);
        }
    }
    tgaFlipVertically(image);    
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

void swapf(double *a, double *b){
    double t = *a;
    *a = *b;
    *b = t;    
}

int abs(int a) {
    return (a >= 0) ? a : -a;
}

int sign(int a) {
    return (a < 0) ? -1 : 1;
}
