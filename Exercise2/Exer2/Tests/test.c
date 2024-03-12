 
#include <stdio.h> 
#include <stdlib.h> 
#include <pthread.h> 
  
// A normal C function that is executed as a thread  
// when its name is specified in pthread_create() 
void *myThreadFun(void *vargp) 
{ 
    // printf("Printing GeeksQuiz from Thread \n"); 
    int *ptr = (int *)vargp;
    printf("Thread ID: %d\n", *ptr);
    return NULL; 
} 
   
int main() 
{ 
    pthread_t thread_id; 
    printf("Before Thread\n"); 
    for(int i=0; i<5; i++){
        int *newID = malloc(sizeof(int));
        *newID = i;
        pthread_create(&thread_id, NULL, myThreadFun, newID);
    }
    pthread_join(thread_id, NULL);
    
    pthread_join(thread_id, NULL); 
    printf("After Thread\n"); 
    exit(0); 
}