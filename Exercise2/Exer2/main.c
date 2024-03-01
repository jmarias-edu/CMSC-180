#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <pthread.h> 
#include <sys/time.h>

struct pearson_in {
    int n;
    double* y;
    double* matrix;
} PEAR_IN;

struct pearson_out {
    double* output;
} PEAR_OUT;

struct y_info {
    double ybar;
    double ysquaresum;
} y_info;

struct thread_args {
    int id;
    int start;
    int end;
    int extra;
} targs;

struct x_info {
    double *xbar;
    double *xysum;
    double *xsquaresum;
} x_info;

struct pearson_in *tocalc;
struct pearson_out *output;
struct y_info *yinfo;
struct x_info *xinfo;

void populate(struct pearson_in *in, struct pearson_out *out, int x){
    in->n = x;
    in->y = (double *)(((char *)in) + sizeof(struct pearson_in));
    in->matrix = (double *)(((char *)in) + sizeof(struct pearson_in) + x*sizeof(double));
    out->output = (double *)(((char *)out) + sizeof(struct pearson_out));

    // double weight[] = {3.63, 3.02, 3.82, 3.42, 3.59, 2.87, 3.03, 3.46, 3.36, 3.3};
    // double length[] = {53.1, 49.7, 48.4, 54.2, 54.9, 43.7, 47.2, 45.2, 54.4, 50.4};

    for(int i=0; i<x; i++){
        in->y[i] = (double) (rand() % 10000 + 1 ) / 100;
        // in->y[i] = length[i];
    }

    for(int i=0; i<x; i++){
        for(int j=0; j<x; j++){
            in->matrix[i*x+j] = (double) (rand() % 10000 + 1 ) / 100;
            // in->matrix[i*x+j] = weight[i];
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

void pearsonCor(){

    double xbar=0, ybar=0, xysum=0, xsquaresum=0, ysquaresum=0;

    for(int a=0; a<tocalc->n; a++){
        ybar += tocalc->y[a];
    }
    ybar = ybar/tocalc->n;

    for(int a=0; a<tocalc->n; a++){
        ysquaresum += pow(tocalc->y[a]-ybar, 2);
    }

    for(int j=0; j<tocalc->n; j++){
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for(int a=0; a<tocalc->n; a++){
            xbar += tocalc->matrix[a*tocalc->n + j];
        }
        xbar = xbar/tocalc->n;

        for(int a=0; a<tocalc->n; a++){
            xysum = xysum + (tocalc->matrix[a*tocalc->n + j] - xbar)*(tocalc->y[a] - ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a*tocalc->n + j] - xbar), 2);
        }

        output->output[j] = xysum / pow((xsquaresum*ysquaresum), 0.5);
    }
}

void *pearsonCorThread(void *vargp){
    struct thread_args *args = (struct thread_args *)vargp;
    printf("Thread ID: %d, Residuals: %d\n", args->id, args->extra);
    
    double xbar=0, xysum=0, xsquaresum=0;

    for(int j=args->start; j<args->end; j++){
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for(int a=0; a<tocalc->n; a++){
            xbar += tocalc->matrix[a*tocalc->n + j];
        }
        xbar = xbar/tocalc->n;

        for(int a=0; a<tocalc->n; a++){
            xysum = xysum + (tocalc->matrix[a*tocalc->n + j] - xbar)*(tocalc->y[a] - yinfo->ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a*tocalc->n + j] - xbar), 2);
        }

        output->output[j] = xysum / pow((xsquaresum*yinfo->ysquaresum), 0.5);

    }
    if(args->extra>0){
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for(int a=0; a<tocalc->n; a++){
            xbar += tocalc->matrix[a*tocalc->n + args->extra];
        }
        xbar = xbar/tocalc->n;

        for(int a=0; a<tocalc->n; a++){
            xysum = xysum + (tocalc->matrix[a*tocalc->n + args->extra] - xbar)*(tocalc->y[a] - yinfo->ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a*tocalc->n + args->extra] - xbar), 2);
        }

        output->output[args->extra] = xysum / pow((xsquaresum*yinfo->ysquaresum), 0.5);
    }
    

    printf("Thread ID %d finished\n", args->id);
    return NULL;
}

void pearsonCorMult(int t){
    pthread_t *tid=malloc(sizeof(pthread_t)*t);
    struct thread_args **threadinfos = malloc(sizeof(struct thread_args *)*t);

    int size=tocalc->n/t;
    int hasResidual = tocalc->n%t!=0;
    
    yinfo = malloc(sizeof(struct y_info));

    for(int a=0; a<tocalc->n; a++){
        yinfo->ybar += tocalc->y[a];
    }
    yinfo->ybar = yinfo->ybar/tocalc->n;

    for(int a=0; a<tocalc->n; a++){
        yinfo->ysquaresum += pow(tocalc->y[a]-yinfo->ybar, 2);
    }

    for(int i=0; i<t; i++){
        struct thread_args *newThreadInfo = malloc(sizeof(struct thread_args));
        threadinfos[i] = newThreadInfo;
        newThreadInfo->id = i;
        newThreadInfo->start = i*size;
        newThreadInfo->end = (i+1)*size;
        if(hasResidual && (size*t+i<tocalc->n)){
            newThreadInfo->extra = size*t+i;
        } else {
            newThreadInfo->extra = -1;
        }
        pthread_create(&tid[i], NULL, pearsonCorThread, (void *)newThreadInfo);
        printf("Thread ID: %d created\n", newThreadInfo->id);
    }
    for(int i=0; i<t; i++){
        pthread_join(tid[i], NULL);
    }
    free(yinfo);
    for(int x=0; x<t; x++){
        free(threadinfos[x]);
    }
    free(threadinfos);
    return;
}

void *pearsonCorThreadRow(void *vargp){
    struct thread_args *args = (struct thread_args *)vargp;
    printf("Thread ID: %d, Residuals: %d\n", args->id, args->extra);
    
    double xbar=0, xysum=0, xsquaresum=0;

    for(int j=args->start; j<args->end; j++){
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for(int a=0; a<tocalc->n; a++){
            xbar += tocalc->matrix[a*tocalc->n + j];
        }
        xbar = xbar/tocalc->n;

        for(int a=0; a<tocalc->n; a++){
            xysum = xysum + (tocalc->matrix[a*tocalc->n + j] - xbar)*(tocalc->y[a] - yinfo->ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a*tocalc->n + j] - xbar), 2);
        }

        output->output[j] = xysum / pow((xsquaresum*yinfo->ysquaresum), 0.5);

    }
    if(args->extra>0){
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for(int a=0; a<tocalc->n; a++){
            xbar += tocalc->matrix[a*tocalc->n + args->extra];
        }
        xbar = xbar/tocalc->n;

        for(int a=0; a<tocalc->n; a++){
            xysum = xysum + (tocalc->matrix[a*tocalc->n + args->extra] - xbar)*(tocalc->y[a] - yinfo->ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a*tocalc->n + args->extra] - xbar), 2);
        }

        
    }
    

    printf("Thread ID %d finished\n", args->id);
    return NULL;
}


void pearsonCorMultRow(int t){
    pthread_t *tid=malloc(sizeof(pthread_t)*t);
    struct thread_args **threadinfos = malloc(sizeof(struct thread_args *)*t);

    int size=tocalc->n/t;
    int hasResidual = tocalc->n%t!=0;
    
    yinfo = malloc(sizeof(y_info));

    for(int a=0; a<tocalc->n; a++){
        yinfo->ybar += tocalc->y[a];
    }
    yinfo->ybar = yinfo->ybar/tocalc->n;

    for(int a=0; a<tocalc->n; a++){
        yinfo->ysquaresum += pow(tocalc->y[a]-yinfo->ybar, 2);
    }

    for(int i=0; i<t; i++){
        struct thread_args *newThreadInfo = malloc(sizeof(struct thread_args));
        threadinfos[i] = newThreadInfo;
        newThreadInfo->id = i;
        newThreadInfo->start = i*size;
        newThreadInfo->end = (i+1)*size;
        if(hasResidual && (size*t+i<tocalc->n)){
            newThreadInfo->extra = size*t+i;
        } else {
            newThreadInfo->extra = -1;
        }
        pthread_create(&tid[i], NULL, pearsonCorThread, (void *)newThreadInfo);
        printf("Thread ID: %d created\n", newThreadInfo->id);
    }
    for(int i=0; i<t; i++){
        pthread_join(tid[i], NULL);
    }

    for(int i=0; i<tocalc; i++){
        output->output[i] = xinfo->xysum[i]  / pow((xinfo->xsquaresum[i]*yinfo->ysquaresum), 0.5);
    }
        

    free(yinfo);
    for(int x=0; x<t; x++){
        free(threadinfos[x]);
    }
    free(threadinfos);
    return;
}

int main(){
    
    // Initializing clock variables
    double t;
    struct timeval start, end;
    
    // Initialize Randomizer
    time_t r;
    srand((unsigned) time(&r));

    // Initialize N and tn
    int n = 0;
    int tn = 0;

    // Initialize File Pointer
    FILE *fptr;

    while(1){
        scanf("%d %d", &n, &tn);
        fptr = fopen("output.txt", "a");
        
        if(n<0 || tn<0){break;} //exits when n or tn is negative

        tocalc = malloc(sizeof(PEAR_IN)+ n*n*sizeof(double) + n*sizeof(double));
        output = malloc(sizeof(PEAR_OUT) + n*sizeof(double));
        populate(tocalc, output, n);
        // printInput(tocalc);

        // Start Tracking time
        gettimeofday(&start, NULL);
        printf("Clock Started\n");
        // pearsonCor(tocalc, output);
        pearsonCorMult(tn);
        // printOutput(output, n);

        // Checking for final runtime
        gettimeofday(&end, NULL);
        t = (end.tv_sec - start.tv_sec) * 1e6;
        t = (t + (end.tv_usec - start.tv_usec)) * 1e-6;

        // Printing to output file and terminal
        printf("n: %d # of Threads: %d Time elapsed: %0.5f\n", n, tn, t);
        fprintf(fptr, "n: %d # of Threads: %d Time elapsed: %0.5f\n", n, tn, t);

        // Free up the memory
        free(tocalc);
        free(output);
        fclose(fptr);
    }
}