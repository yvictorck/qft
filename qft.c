// Author: Wes Kendall
// Copyright 2012 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header in tact.
//
// Program that computes the average of an array of elements in parallel using
// MPI_Scatter and MPI_Gather
//
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
  
// Creates an array of random numbers. Each number has a value from 0 - 1
float *create_rand_nums(int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)RAND_MAX);

  }
  return rand_nums;
}


float *create_nums(int num_elements,int x) {
  float *nums = (float *)malloc(sizeof(float) * num_elements);
  assert(nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    nums[i] = x;
        printf("create num:%f\n",nums[i] );
  }
  return nums;
}
// Computes the average of an array of numbers
float compute_avg(float *array, int num_elements) {
  float sum = 0.f;
  int i;
  for (i = 0; i < num_elements; i++) {
    sum += array[i];
  }
  return sum / num_elements;
}
// Computes the average of an array of numbers
float *compute_neg(float *array, int num_elements,int n) {
  float *nums = (float *)malloc(sizeof(float) * num_elements);
  assert(nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    // printf("%f\n",array[i] );
    nums[i] = array[i]+n;
    printf("world_rank: %d with number %f\n",n,nums[i] );
  }
  return nums;
}


int main(int argc, char** argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: avg num_elements_per_proc\n");
    exit(1);
  }

  int num_elements_per_proc = atoi(argv[1]);
  // Seed the random number generator to get different results each time
  srand(time(NULL));

  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  printf("world_size: %d\n",world_size );

  // Create a random array of elements on the root process. Its total
  // size will be the number of elements per process times the number
  // of processes
  float *rand_nums = NULL;
  if (world_rank == 0) {
    rand_nums = create_nums(num_elements_per_proc * world_size,world_rank);
  }

  // For each process, create a buffer that will hold a subset of the entire
  // array
  float *sub_rand_nums = (float *)malloc(sizeof(float) * num_elements_per_proc);
  assert(sub_rand_nums != NULL);

  // Scatter the random numbers from the root process to all processes in
  // the MPI world
  MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums,
              num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Compute the average of your subset
  float *rec_indi_nums = NULL;
  rec_indi_nums = compute_neg(sub_rand_nums, num_elements_per_proc,world_rank);
 
  // Gather all partial averages down to the root process
  float *rec_gather_nums = NULL;
  if (world_rank == 0) {
    rec_gather_nums = (float *)malloc(sizeof(float) * num_elements_per_proc* world_size);
    assert(rec_gather_nums != NULL);
  }
  MPI_Gather(rec_indi_nums, num_elements_per_proc, MPI_FLOAT, rec_gather_nums, num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Now that we have all of the partial averages on the root, compute the
  // total average of all numbers. Since we are assuming each process computed
  // an average across an equal amount of elements, this computation will
  // produce the correct answer.
  if (world_rank == 0) {
    //float avg = compute_avg(sub_avgs, world_size);
    // printf("Avg of all elements is %f\n", avg);
    // // Compute the average across the original data for comparison
    // float original_data_avg =
    //   compute_avg(rand_nums, num_elements_per_proc * world_size);
    // printf("Avg computed across original data is %f\n", original_data_avg);
    int i;
  for (i = 0; i < num_elements_per_proc*world_size; i++) {
    printf("%f ", rec_gather_nums[i]);
      printf( "\n" );
  }



  }

  // Clean up
  if (world_rank == 0) {
    free(rec_indi_nums);
    free(rec_gather_nums);
  }
  free(sub_rand_nums);
 
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
