#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tga.h"
#include "model.h"

void swap(int *a, int *b);
int abs(int a);
int sign(int a);

// Triangle rasterization algorithm using a scan line
void triangle(tgaImage *image,
              int x0, int y0,
              int x1, int y1,
              int x2, int y2,
              tgaColor color);

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

        image = tgaNewImage(4000, 4000, RGB);
        if (!image) {
            perror("tgaNewImage");
            rv = -1;
            break;
        }
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
              int x0, int y0,
              int x1, int y1,
              int x2, int y2,
              tgaColor color)
{
    // Sorting vertices by 'y' coord
    if (y0 > y1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }
    if (y0 > y2) {
        swap(&x0, &x2);
        swap(&y0, &y2);
    }
    if (y1 > y2) {
        swap(&x1, &x2);
        swap(&y1, &y2);
    }

    // Computing x coord deltas
    double dx01 = 0, dx02 = 0, dx12 = 0;
    if (y0 != y1) {
        dx01 = x1 - x0;
        dx01 /= y1 - y0;
    }
    if (y0 != y2) {
        dx02 = x2 - x0;
        dx02 /= y2 - y0;
    }
    if (y1 != y2) {
        dx12 = x2 - x1;
        dx12 /= y2 - y1;
    }
    double _dx02 = dx02;

    if (dx01 > dx02) {
        double t = dx01;
        dx01 = dx02;
        dx02 = t;
    }

    if (dx12 > _dx02) {
        double t = dx12;
        dx12 = _dx02;
        _dx02 = t;
    }

    //Fill the triangle
    int y, x;
    double xleft = x0, xright = x0;
    for (y = y0; y <= y2; ++y) {
        for (x = xleft; x <= xright; ++x) {
            tgaSetPixel(image, x, y, color);
        }
        xleft += (y < y1) ? dx01 : _dx02;
        xright += (y < y1) ? dx02 : dx12;
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
    int light[3] = {0, 1, 0};

    for (face = 0; face < model->nface; ++face) {

        // Translating into screen coordinates
        for (vert = 0; vert < 3; ++vert) {
            p[vert] = getVertex(model, face, vert);
            screen_coords[vert][0] = ((*p[vert])[0] / 4000 + 1.0) * w / 2;
            screen_coords[vert][1] = ((*p[vert])[2] / 4000 + 1.0) * h / 2;
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
        if(intensity < 0){
            // Fill the triangle
            triangle(image, screen_coords[0][0], screen_coords[0][1],
                        screen_coords[1][0], screen_coords[1][1],
                        screen_coords[2][0], screen_coords[2][1],
                        tgaRGB((int)(fabs(intensity) * 255) % 255, (int)(fabs(intensity) * 255) % 255, (int)(fabs(intensity) * 255) % 255));
        }
    }
    tgaFlipVertically(image);    
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int abs(int a) {
    return (a >= 0) ? a : -a;
}

int sign(int a) {
    return (a < 0) ? -1 : 1;
}
