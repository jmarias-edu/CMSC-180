#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

struct pearson_in {
    int n;
    double* y;
    double* matrix;
} PEAR_IN;

struct pearson_out {
    double* output;
} PEAR_OUT;

void populate(struct pearson_in *in, struct pearson_out *out, int x){
    in->n = x;
    in->y = (double *)(((char *)in) + sizeof(struct pearson_in));
    in->matrix = (double *)(((char *)in) + sizeof(struct pearson_in) + x*sizeof(double));
    out->output = (double *)(((char *)out) + sizeof(struct pearson_out));

    double weight[] = {3.63, 3.02, 3.82, 3.42, 3.59, 2.87, 3.03, 3.46, 3.36, 3.3};
    double length[] = {53.1, 49.7, 48.4, 54.2, 54.9, 43.7, 47.2, 45.2, 54.4, 50.4};

    for(int i=0; i<x; i++){
        // in->y[i] = (double) (rand() % 10000 + 1 ) / 100;
        in->y[i] = length[i];
    }

    for(int i=0; i<x; i++){
        for(int j=0; j<x; j++){
            // in->matrix[i*x+j] = (double) (rand() % 10000 + 1 ) / 100;
            in->matrix[i*x+j] = weight[i];
        }
    }
}

void printInput(struct pearson_in *in){
    printf("n: %f\n", in->n);

    printf("y array:\n");
    for(int i=0; i<in->n; i++){
        printf("%f ", in->y[i]);
    }
    printf("\n");

    printf("x matrix:\n");
    for(int i=0; i<in->n; i++){
        for(int j=0; j<in->n; j++){
            printf("%f ", in->matrix[i*in->n+j]);
        }
        printf("\n");
    }
}

void printOutput(struct pearson_out *out, int n){
    printf("R values:\n");
    for(int i=0; i<n; i++){
        printf("%f ", out->output[i]);
    }
    printf("\n");
}

void pearson_cor(struct pearson_in *in, struct pearson_out *out){

    double xbar=0, ybar=0, xysum=0, xsquaresum=0, ysquaresum=0;

    for(int a=0; a<in->n; a++){
        ybar += in->y[a];
    }
    ybar = ybar/in->n;

    for(int a=0; a<in->n; a++){
        ysquaresum += pow(in->y[a]-ybar, 2);
    }

    for(int j=0; j<in->n; j++){
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for(int a=0; a<in->n; a++){
            xbar += in->matrix[a*in->n + j];
        }
        xbar = xbar/in->n;

        for(int a=0; a<in->n; a++){
            xysum = xysum + (in->matrix[a*in->n + j] - xbar)*(in->y[a] - ybar);
            xsquaresum = xsquaresum + pow((in->matrix[a*in->n + j] - xbar), 2);
        }

        out->output[j] = xysum / pow((xsquaresum*ysquaresum), 0.5);
    }
}

int main(){
    
    // Initializing clock variables
    double t;
    struct timeval start, end;
    
    // Initialize Randomizer
    time_t r;
    srand((unsigned) time(&r));

    // Initialize N
    int n = 0;

    // Allocate Memory for Pearson
    struct pearson_in *tocalc = malloc(sizeof(PEAR_IN)+ n*n*sizeof(double) + n*sizeof(double));
    struct pearson_out *output = malloc(sizeof(PEAR_OUT) + n*sizeof(double));

    // Initialize File Pointer
    FILE *fptr;

    while(1){
        scanf("%d", &n);
        fptr = fopen("output.txt", "a");
        
        if(n<0){break;}
        struct pearson_in *tocalc = malloc(sizeof(PEAR_IN)+ n*n*sizeof(double) + n*sizeof(double));
        struct pearson_out *output = malloc(sizeof(PEAR_OUT) + n*sizeof(double));
        populate(tocalc, output, n);
        // printInput(tocalc);

        printf("Clock Started\n");
        gettimeofday(&start, NULL);
        pearson_cor(tocalc, output);

        printOutput(output, n);

        // Checking for final runtime

        gettimeofday(&end, NULL);
        t = (end.tv_sec - start.tv_sec) * 1e6;
        t = (t + (end.tv_usec - start.tv_usec)) * 1e-6;
        printf("n: %d Time elapsed: %0.5f\n", n, t);
        fprintf(fptr, "n: %d Time elapsed: %0.5f\n", n, t);
        free(tocalc);
        free(output);
        fclose(fptr);
    }
}

