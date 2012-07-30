#include "iio.c"
#include <stdio.h> // only for "fprintf"
#include <stdlib.h> // only for "free"

int main(int c, char *v[]){
    
    if (c!=3){
        printf("Usage: %s input_image output_image\n",v[0]);
    } else {
        int w, h, pixeldim;
        float *x = iio_read_image_float_vec(v[1], &w, &h, &pixeldim);
        fprintf(stderr, "Got a %dx%d image with %d channels\n", w, h, pixeldim);
        
        // some processing here
        float *xgray = malloc(w*h*sizeof(float));
        if (pixeldim==3){
                int i;
                for(i=0;i<w*h;i++){
                        xgray[i] =  (6968*x[pixeldim*i] + 23434*x[pixeldim*i + 1] + 2366*x[pixeldim*i + 2])/32768;
                }
        // Y = (6968 R + 23434 G + 2366 B) / 32768
        }
        else {
            xgray = x;
        }

        printf("w = %d, h = %d, pixeldim = %d.\n",w,h,pixeldim);
        
        iio_save_image_float_vec(v[2], xgray, w, h, 1);
        free(x);
        return 0;
    }
    
}
