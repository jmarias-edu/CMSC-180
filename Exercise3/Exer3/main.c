// Jarem Thimoty M. Arias
// 2020-00782
// CMSC 180 - T5L
// Laboratory Problem 3 Code

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include <sched.h>

//structure for pearson input
struct pearson_in
{
    int n;
    double *y;
    double *matrix;
} PEAR_IN;

//structure for pearson output
struct pearson_out
{
    double *output;
} PEAR_OUT;

//structure for y array related information
struct y_info
{
    double ybar;
    double ysquaresum;
} y_info;

//structure to pass to threads as arguments for range over matrix
struct thread_args
{
    int id;
    int start;
    int end;
    int extra;
} targs;

//structure for x related information (used in row submatrice)
struct x_info
{
    double *xbar;
    int *xbarcount;
    int xbarflag;
    double *xysum;
    double *xsquaresum;
} x_info;

//initiating global structures
struct pearson_in *tocalc;
struct pearson_out *output;
struct y_info *yinfo;
struct x_info *xinfo;
pthread_mutex_t *locklist;

//populates x matrix
void populate(struct pearson_in *in, struct pearson_out *out, int x)
{
    in->n = x;
    in->y = (double *)(((char *)in) + sizeof(struct pearson_in));
    in->matrix = (double *)(((char *)in) + sizeof(struct pearson_in) + x * sizeof(double));
    out->output = (double *)(((char *)out) + sizeof(struct pearson_out));

    for (int i = 0; i < x; i++)
    {
        in->y[i] = (double)(rand() % 10000 + 1) / 100;
    }

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < x; j++)
        {
            in->matrix[i * x + j] = (double)(rand() % 10000 + 1) / 100;
        }
    }
}

//populates matrix with known values (From handout)
void testPopulate(struct pearson_in *in, struct pearson_out *out, int x)
{
    in->n = x;
    in->y = (double *)(((char *)in) + sizeof(struct pearson_in));
    in->matrix = (double *)(((char *)in) + sizeof(struct pearson_in) + x * sizeof(double));
    out->output = (double *)(((char *)out) + sizeof(struct pearson_out));

    double weight[] = {3.63, 3.02, 3.82, 3.42, 3.59, 2.87, 3.03, 3.46, 3.36, 3.3};
    double length[] = {53.1, 49.7, 48.4, 54.2, 54.9, 43.7, 47.2, 45.2, 54.4, 50.4};

    for (int i = 0; i < x; i++)
    {
        in->y[i] = length[i];
    }

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < x; j++)
        {
            in->matrix[i * x + j] = weight[i];
        }
    }
}

//accessory function to print input
void printInput(struct pearson_in *in)
{
    printf("n: %f\n", in->n);

    printf("y array:\n");
    for (int i = 0; i < in->n; i++)
    {
        printf("%f ", in->y[i]);
    }
    printf("\n");

    printf("x matrix:\n");
    for (int i = 0; i < in->n; i++)
    {
        for (int j = 0; j < in->n; j++)
        {
            printf("%f ", in->matrix[i * in->n + j]);
        }
        printf("\n");
    }
}

//accessory function to print output
void printOutput(struct pearson_out *out, int n)
{
    printf("R values:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%f ", out->output[i]);
    }
    printf("\n");
}

//serial pearson correlation
void pearsonCor()
{
    //initiates variables
    double xbar = 0, ybar = 0, xysum = 0, xsquaresum = 0, ysquaresum = 0;

    //calculates for average of y
    for (int a = 0; a < tocalc->n; a++)
    {
        ybar += tocalc->y[a];
    }
    ybar = ybar / tocalc->n;

    //calculates for summation of squares
    for (int a = 0; a < tocalc->n; a++)
    {
        ysquaresum += pow(tocalc->y[a] - ybar, 2);
    }

    //main computation
    for (int j = 0; j < tocalc->n; j++)
    {
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        //computes for average of x
        for (int a = 0; a < tocalc->n; a++)
        {
            xbar += tocalc->matrix[a * tocalc->n + j];
        }
        xbar = xbar / tocalc->n;

        //computes for summation of x and y along with summation of x squared
        for (int a = 0; a < tocalc->n; a++)
        {
            xysum = xysum + (tocalc->matrix[a * tocalc->n + j] - xbar) * (tocalc->y[a] - ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a * tocalc->n + j] - xbar), 2);
        }

        //computes the final output
        output->output[j] = xysum / pow((xsquaresum * ysquaresum), 0.5);
    }
}

// function of by column submatrix thread
void *pearsonCorThread(void *vargp)
{
    //Gets thread arguments
    struct thread_args *args = (struct thread_args *)vargp;
    printf("Thread ID: %d, Residuals: %d\n", args->id, args->extra);

    double xbar = 0, xysum = 0, xsquaresum = 0;

    //same computation as serial above, only difference is a determined range
    for (int j = args->start; j < args->end; j++)
    {
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for (int a = 0; a < tocalc->n; a++)
        {
            xbar += tocalc->matrix[a * tocalc->n + j];
        }
        xbar = xbar / tocalc->n;

        for (int a = 0; a < tocalc->n; a++)
        {
            xysum = xysum + (tocalc->matrix[a * tocalc->n + j] - xbar) * (tocalc->y[a] - yinfo->ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a * tocalc->n + j] - xbar), 2);
        }

        output->output[j] = xysum / pow((xsquaresum * yinfo->ysquaresum), 0.5);
    }
    //same as above, only computes for assigned extra column
    if (args->extra > 0)
    {
        xbar = 0;
        xysum = 0;
        xsquaresum = 0;

        for (int a = 0; a < tocalc->n; a++)
        {
            xbar += tocalc->matrix[a * tocalc->n + args->extra];
        }
        xbar = xbar / tocalc->n;

        for (int a = 0; a < tocalc->n; a++)
        {
            xysum = xysum + (tocalc->matrix[a * tocalc->n + args->extra] - xbar) * (tocalc->y[a] - yinfo->ybar);
            xsquaresum = xsquaresum + pow((tocalc->matrix[a * tocalc->n + args->extra] - xbar), 2);
        }

        output->output[args->extra] = xysum / pow((xsquaresum * yinfo->ysquaresum), 0.5);
    }

    printf("Thread ID %d finished\n", args->id);
    return NULL;
}

// main function for column submatrix thread
void pearsonCorMult(int t, int affinity, int cores)
{
    //initiates thread variables
    pthread_t *tid = malloc(sizeof(pthread_t) * t);
    struct thread_args **threadinfos = malloc(sizeof(struct thread_args *) * t);
    cpu_set_t mask;

    //computes for distribution of submatrices through threads
    int size = tocalc->n / t;
    int hasResidual = tocalc->n % t != 0;

    //precomputes for info regarding y to be used by all threads
    yinfo = malloc(sizeof(struct y_info));

    for (int a = 0; a < tocalc->n; a++)
    {
        yinfo->ybar += tocalc->y[a];
    }
    yinfo->ybar = yinfo->ybar / tocalc->n;

    for (int a = 0; a < tocalc->n; a++)
    {
        yinfo->ysquaresum += pow(tocalc->y[a] - yinfo->ybar, 2);
    }

    //initiates threads
    for (int i = 0; i < t; i++)
    {
        //creates new threads and assigns data
        struct thread_args *newThreadInfo = malloc(sizeof(struct thread_args));
        threadinfos[i] = newThreadInfo;
        newThreadInfo->id = i;
        newThreadInfo->start = i * size;
        newThreadInfo->end = (i + 1) * size;
        if (hasResidual && (size * t + i < tocalc->n)) //adds extra column
        {
            newThreadInfo->extra = size * t + i;
        }
        else
        {
            newThreadInfo->extra = -1;
        }
        if(affinity==1) // if core affinity is activated, assign core to thread
        {
            CPU_ZERO(&mask);
            CPU_SET(i%cores, &mask);
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &mask);
            pthread_create(&tid[i], &attr, pearsonCorThread, (void *)newThreadInfo);
        }
        else // if core affinity is not activated
        {
            pthread_create(&tid[i], NULL, pearsonCorThread, (void *)newThreadInfo);
        }
        printf("Thread ID: %d created\n", newThreadInfo->id);
    }
    for (int i = 0; i < t; i++) //joins the threads
    {
        pthread_join(tid[i], NULL);
    }
    free(yinfo); //frees all the information
    for (int x = 0; x < t; x++)
    {
        free(threadinfos[x]);
    }
    free(threadinfos);
    return;
}

//function for row submatrix threads
void *pearsonCorThreadRow(void *vargp)
{   //gets thread arguments
    struct thread_args *args = (struct thread_args *)vargp;
    printf("Thread ID: %d, Residuals: %d\n", args->id, args->extra);

    //adds to xbar array in xinfo for computing xbar
    for (int j = args->start; j < args->end; j++)
    {
        for (int a = 0; a < tocalc->n; a++)
        {
            pthread_mutex_lock(&locklist[a]);
            xinfo->xbar[a] += tocalc->matrix[j * tocalc->n + a];
            xinfo->xbarcount[a]++;
            pthread_mutex_unlock(&locklist[a]);
        }
    }
    //computes for extra row assigned
    if (args->extra > 0)
    {
        for (int a = 0; a < tocalc->n; a++)
        {
            pthread_mutex_lock(&locklist[a]);
            xinfo->xbar[a] += tocalc->matrix[args->extra * tocalc->n + a];
            xinfo->xbarcount[a]++;
            pthread_mutex_unlock(&locklist[a]);
        }
    }

    //waits for main thread to finalize xbar computation
    printf("Thread ID %d Waiting\n", args->id);
    while (xinfo->xbarflag == 0); // waits for main thread to finish computing xbar
    printf("Thread ID %d Continuing\n", args->id);

    //continues computation that requires xbar value
    for (int j = args->start; j < args->end; j++)
    {
        for (int a = 0; a < tocalc->n; a++)
        {
            pthread_mutex_lock(&locklist[a]);
            xinfo->xysum[a] = xinfo->xysum[a] + (tocalc->matrix[j * tocalc->n + a] - xinfo->xbar[a]) * (tocalc->y[j] - yinfo->ybar);
            xinfo->xsquaresum[a] = xinfo->xsquaresum[a] + pow((tocalc->matrix[j * tocalc->n + a] - xinfo->xbar[a]), 2);
            pthread_mutex_unlock(&locklist[a]);
        }
    }
    //computes for extra row
    if (args->extra > 0)
    {
        for (int a = 0; a < tocalc->n; a++)
        {
            pthread_mutex_lock(&locklist[a]);
            xinfo->xysum[a] = xinfo->xysum[a] + (tocalc->matrix[args->extra * tocalc->n + a] - xinfo->xbar[a]) * (tocalc->y[args->extra] - yinfo->ybar);
            xinfo->xsquaresum[a] = xinfo->xsquaresum[a] + pow((tocalc->matrix[args->extra * tocalc->n + a] - xinfo->xbar[a]), 2);
            pthread_mutex_unlock(&locklist[a]);
        }
    }

    printf("Thread ID %d finished\n", args->id);
    return NULL;
}

//main function for row submatrix
void pearsonCorMultRow(int t, int affinity, int cores)
{   //initializes thread variables
    pthread_t *tid = malloc(sizeof(pthread_t) * t);
    struct thread_args **threadinfos = malloc(sizeof(struct thread_args *) * t);
    cpu_set_t mask;

    //computes for distribution of submatrices through threads
    int size = tocalc->n / t;
    int hasResidual = tocalc->n % t != 0;
    //assigns information to xinfo to be shared among threads
    yinfo = malloc(sizeof(y_info));
    xinfo = malloc(sizeof(struct x_info) + tocalc->n * 3 * sizeof(double) + tocalc->n * sizeof(int));
    memset(xinfo, 0, sizeof(struct x_info) + tocalc->n * 3 * sizeof(double) + tocalc->n * sizeof(int));
    xinfo->xbar = (double *)(((char *)xinfo) + sizeof(struct x_info));
    xinfo->xysum = (double *)(((char *)xinfo) + sizeof(struct x_info) + tocalc->n * sizeof(double));
    xinfo->xsquaresum = (double *)(((char *)xinfo) + sizeof(struct x_info) + 2 * tocalc->n * sizeof(double));
    xinfo->xbarcount = (int *)(((char *)xinfo) + sizeof(struct x_info) + 3 * tocalc->n * sizeof(double));
    xinfo->xbarflag = 0;
    locklist = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t) * tocalc->n); //initiates array of mutex

    //populates mutex to be used by threads to prevent race conditions
    for(int i = 0; i<tocalc->n; i++){
        pthread_mutex_init(&locklist[i], NULL);
    }

    //precomputes for y array related info
    for (int a = 0; a < tocalc->n; a++)
    {
        yinfo->ybar += tocalc->y[a];
    }
    yinfo->ybar = yinfo->ybar / tocalc->n;

    for (int a = 0; a < tocalc->n; a++)
    {
        yinfo->ysquaresum += pow(tocalc->y[a] - yinfo->ybar, 2);
    }

    //initiates threads
    for (int i = 0; i < t; i++)
    {
        struct thread_args *newThreadInfo = malloc(sizeof(struct thread_args));
        threadinfos[i] = newThreadInfo;
        newThreadInfo->id = i;
        newThreadInfo->start = i * size;
        newThreadInfo->end = (i + 1) * size;
        if (hasResidual && (size * t + i < tocalc->n))
        {
            newThreadInfo->extra = size * t + i;
        }
        else
        {
            newThreadInfo->extra = -1;
        }
        
        if(affinity==1) // if core affinity is activated, assign core to thread
        {
            CPU_ZERO(&mask);
            CPU_SET(i%cores, &mask);
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &mask);
            pthread_create(&tid[i], &attr, pearsonCorThreadRow, (void *)newThreadInfo);
        }
        else // if core affinity is not activated
        {
            pthread_create(&tid[i], NULL, pearsonCorThreadRow, (void *)newThreadInfo);
        }
        printf("Thread ID: %d created\n", newThreadInfo->id);
    }

    printf("Main Thread waiting for completion of xbar\n");

    //loop to wait for threads to finish adding all x
    restartloop:
    for (int i = 0; i < tocalc->n; i++)
    {
        if (xinfo->xbarcount[i] != tocalc->n)
        {
            goto restartloop;
        }
    }

    //computes for final value of xbar
    printf("Main thread computing for final value of xbar\n");
    for (int i = 0; i < tocalc->n; i++)
    {
        xinfo->xbar[i] = xinfo->xbar[i] / tocalc->n;
    }

    //lets threads continue working after computing xbar
    xinfo->xbarflag = 1;
    printf("xbar computed, main thread will continue waiting\n");

    for (int i = 0; i < t; i++) //joins threads
    {
        pthread_join(tid[i], NULL);
    }

    for (int i = 0; i < tocalc->n; i++) //computes for final values
    {
        output->output[i] = xinfo->xysum[i] / pow((xinfo->xsquaresum[i] * yinfo->ysquaresum), 0.5);
    }

    free(yinfo); //frees information
    free(xinfo);
    free(locklist);
    for (int x = 0; x < t; x++)
    {
        free(threadinfos[x]);
    }
    free(threadinfos);
    return;
}

int main()
{
    // Initializing clock variables
    double t;
    struct timeval start, end;

    // Initialize Randomizer
    time_t r;
    srand((unsigned)time(&r));

    // Initialize N and tn
    int n = 0;
    int tn = 0;
    int type = 0;
    int cores = 0;

    // Initialize File Pointer
    FILE *fptr;

    // Set affinity of Main thread to 0
    // cpu_set_t cpus;
    // CPU_ZERO(&cpus);
    // CPU_SET(0, &cpus);
    // sched_setaffinity(0, sizeof(cpu_set_t), &cpus);

    while (1)
    {
        scanf("%d %d %d %d", &n, &tn, &type, &cores);
        fptr = fopen("output.txt", "a");

        if (n < 0 || tn < 0 || type < 0 || cores < 0)
        {
            break;
        } // exits when n, tn, or type is negative

        tocalc = malloc(sizeof(PEAR_IN) + n * n * sizeof(double) + n * sizeof(double));
        output = malloc(sizeof(PEAR_OUT) + n * sizeof(double));
        populate(tocalc, output, n);
        // testPopulate(tocalc, output, n);
        // printInput(tocalc);

        // Start Tracking time
        gettimeofday(&start, NULL);
        printf("Clock Started\n");


        //choosing type of run
        if (type == 0)
        {
            pearsonCor();
        }
        else if (type == 1)
        {
            pearsonCorMult(tn, 0, cores);
        }
        else if (type == 2)
        {
            pearsonCorMultRow(tn, 0, cores);
        }
        else if (type == 3)
        {
            pearsonCorMult(tn, 1, cores);
        }
        else if (type == 4)
        {
            pearsonCorMultRow(tn, 1, cores);
        }

        // Testing Function to check output
        // printOutput(output, n);

        // Checking for final runtime
        gettimeofday(&end, NULL);
        t = (end.tv_sec - start.tv_sec) * 1e6;
        t = (t + (end.tv_usec - start.tv_usec)) * 1e-6;

        // Printing to output file and terminal
        printf("n: %d\t # of Threads: %d\t Type: %d\t # of cores: %d\t Time elapsed: %0.5f\n", n, tn, type, cores, t);
        fprintf(fptr, "n: %d\t # of Threads: %d\t Type %d\t # of cores: %d\t Time elapsed: %0.5f\n", n, tn, type, cores, t);

        // Free up the memory
        free(tocalc);
        free(output);
        fclose(fptr);
    }
}