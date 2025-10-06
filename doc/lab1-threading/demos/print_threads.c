#include <pthread.h>   // Needed for POSIX threads (pthreads)
#include <stdio.h>     // Standard input/output
#include <stdlib.h>    // For malloc() and free()

// Define a struct to hold arguments for each thread
// This lets us pass more than one value to the thread function
typedef struct {
    int num;           // An integer (e.g., thread number)
    const char* msg;   // A pointer to a string (the message)
} ThreadArgs;

// This is the function each thread will run
// It takes a void* argument (can point to anything)
void* print_info(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg; // Cast the void* back to our struct type
    printf("Thread %d says: %s\n", args->num, args->msg); // Print the info
    free(arg); // Free the memory we allocated for this thread's arguments
    return NULL; // Thread returns nothing
}

int main() {
    pthread_t threads[3]; // Array to hold thread IDs (handles)
    const char* messages[] = {"Hello", "from", "threads!"}; // Messages for each thread
    for (int i = 0; i < 3; i++) {
        ThreadArgs* args = malloc(sizeof(ThreadArgs)); // Allocate memory for arguments
        args->num = i + 1; // Set the thread number (1, 2, 3)
        args->msg = messages[i]; // Set the message for this thread
        pthread_create(&threads[i], NULL, print_info, args); // Start the thread, pass args
    }
    for (int i = 0; i < 3; i++) {
        pthread_join(threads[i], NULL); // Wait for each thread to finish
    }
    return 0;
}
