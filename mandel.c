// Cameron Gibson

// CS 4352-001

#include "bitmap.h"
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

typedef struct {
	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double iterationMax;
	int height_begin;
	int height_end;
} thread_args;

int iteration_to_color(int i, int max);
int iterations_at_point(double x, double y, int max);
void* compute_image(void* a);


void show_help()
{
    printf("Use: mandel [options]\n");
    printf("Where options are:\n");
    printf("-n <thread> The number of threads (default=1)\n");
    printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
    printf("-x <coord>  X coordinate of image center point. (default=0)\n");
    printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
    printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
    printf("-W <pixels> Width of the image in pixels. (default=500)\n");
    printf("-H <pixels> Height of the image in pixels. (default=500)\n");
    printf("-o <file>   Set output file. (default=mandel.bmp)\n");
    printf("-h          Show this help text.\n");
    printf("\nSome examples are:\n");
    printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
    printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
    printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main(int argc, char* argv[])
{
    char c;
    // Set default numberThreads to one
    int numberThreads = 1;
    
    // These are the default configuration values used
    // if no command line arguments are given.

    const char* outfile = "mandel.bmp";
    double xcenter = 0;
    double ycenter = 0;
    double scale = 4;
    int image_width = 500;
    int image_height = 500;
    int max = 1000;

    // For each command line argument given,
    // override the appropriate configuration value.

    while ((c = getopt(argc, argv, "x:y:s:W:H:m:n:o:h")) != -1) {
        switch (c) {
        case 'x':
            xcenter = atof(optarg);
            break;
        case 'y':
            ycenter = atof(optarg);
            break;
        case 's':
            scale = atof(optarg);
            break;
        case 'W':
            image_width = atoi(optarg);
            break;
        case 'H':
            image_height = atoi(optarg);
            break;
        case 'm':
            max = atoi(optarg);
            break;
        case 'n':
            numberThreads = atoi(optarg);
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'h':
            show_help();
            exit(1);
            break;
        }
    }

    // Display the configuration of the image.
    printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s numberthreads=%d\n", xcenter, ycenter, scale, max, outfile, numberThreads);
    
    
    // Create an array of threads
    pthread_t* threads = malloc(numberThreads * sizeof(pthread_t));
    
    // A structure for arguments
    thread_args* args = malloc(numberThreads * sizeof(thread_args));

    // Create a bitmap of the appropriate size.
    struct bitmap* bm = bitmap_create(image_width, image_height);

    // Fill it with a dark blue, for debugging
    //bitmap_reset(bm, MAKE_RGBA(0, 0, 255, 0));
    
    int i;
    
    // Divide height into threads
    int heightStep = image_height / numberThreads;
    
    // Allocate a `thread_args`structure and pass 
    // it as the fourth argument to `pthread_create`, which
    // will turn around and pass it to `compute_image`
    
    for(i = 0; i < numberThreads; i++){
    	args[i].bm = bm;
    	args[i].xmin = xcenter - scale;
    	args[i].xmax = xcenter + scale;
    	args[i].ymin = ycenter - scale;
    	args[i].ymax = ycenter + scale;
    	args[i].iterationMax = max;
    	
    	if (i == 0) {
    		args[i].height_begin = 0;
    		args[i].height_end = heightStep;
    	}
    	
    	else{
    		args[i].height_begin = args[i-1].height_end;
    		args[i].height_end = args[i-1].height_end + heightStep;
    	}
    	printf("Creating thread %d\n", i);
    	pthread_create(&threads[i], NULL, compute_image, (void*) &args[i]);
    }
    
    for(i = 0; i < numberThreads; i++){
    	printf("Joining thread %d\n", i);
    	pthread_join(threads[i], NULL);
    }

    // Compute the Mandelbrot image
    //compute_image(bm, xcenter - scale, xcenter + scale, ycenter - scale, ycenter + scale, max);

    // Save the image in the stated file.
    if (!bitmap_save(bm, outfile)) {
        fprintf(stderr, "mandel: couldn't write to %s: %s\n", outfile, strerror(errno));
        return 1;
    }

    return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void* compute_image(void* a)
{
    // `thread_args` contains everything you pass to `compute_image`
    thread_args* args = (thread_args*) a;
    int i, j;

    int width = bitmap_width(args -> bm);
    int height = bitmap_height(args -> bm);
    int start = args -> height_begin;
    int end = args -> height_end;

    // For every pixel in the image...
    for (j = start; j < end; j++) {

        for (i = 0; i < width; i++) {

            // Determine the point in x,y space for that pixel.
            double x = args -> xmin + i * (args -> xmax - args -> xmin) / width;
            double y = args -> ymin + j * (args -> ymax - args -> ymin) / height;

            // Compute the iterations at that point.
            int iters = iterations_at_point(x, y, args -> iterationMax);

            // Set the pixel in the bitmap.
            bitmap_set(args -> bm, i, j, iters);
        }
    }
    return 0;
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point(double x, double y, int max)
{
    double x0 = x;
    double y0 = y;

    int iter = 0;

    while ((x * x + y * y <= 4) && iter < max) {

        double xt = x * x - y * y + x0;
        double yt = 2 * x * y + y0;

        x = xt;
        y = yt;

        iter++;
    }

    return iteration_to_color(iter, max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color(int i, int max)
{
    int gray = 255 * i / max;
    return MAKE_RGBA(gray, gray, gray, 0);
}