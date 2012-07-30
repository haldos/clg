#include "iio.c"
#include <stdio.h> // only for "fprintf"
#include <stdlib.h> // only for "free"



int main(int c, char *v[]){
    
    if (c!=2){
        printf("Usage: %s tiff_image\n",v[0]);
    } else {
        int w, h, pixeldim;
        float *x = iio_read_image_float_vec(v[1], &w, &h, &pixeldim);
        fprintf(stderr, "Got a %dx%d image with %d channels\n", w, h, pixeldim);
    
    float *x2 = malloc(w * h * 2 * sizeof(float));
    int i;
    for(i=0;i<w*h;i++){
        x2[2*i] = -x[2*i + 1];
        x2[2*i + 1] = -x[2*i];
    }
    
    iio_save_image_float_vec(v[1], x2, w, h, 2);
    free(x2);
    }
    
    return 0;
}
