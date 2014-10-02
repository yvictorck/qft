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
// #include <math.h>
#include <cmath>

#include <complex>  
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <bitset>
#include <typeinfo>

double pi;
std::complex<double> I;
namespace ublas = boost::numeric::ublas;

  
// Creates an array of random numbers. Each number has a value from 0 - 1
double *create_rand_nums(int num_elements) {
  double *rand_nums = (double *)malloc(sizeof(double) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (double)RAND_MAX);
  }
  return rand_nums;
}


double *create_nums(int num_elements,int x) {
  double *nums = (double *)malloc(sizeof(double) * num_elements);
  assert(nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    nums[i] = x;
        printf("create num:%f\n",nums[i] );
  }
  return nums;
}

double *create_onelinephases(double arr[],int size) {
  double *nums = (double *)malloc(sizeof(double) * size);
  assert(nums != NULL);
  int i;
  for (i = 0; i < size; i++) {
    nums[i] = arr[i];
  }
  return nums;
}

// len of array = num_point_qft
double *compute_qft(double *array, int n_qft_point, int n) {
  int len_nums_complex = pow(2,n_qft_point);
  std::complex<double> *nums_complex = (std::complex<double> *)malloc(sizeof(std::complex<double>) *len_nums_complex);
  std::complex<double> im(0, 1); 
  
  
  // //
  // for (int i = 0; i <len_nums_complex; i++) 
  // {
  //   nums_complex[i] = std::exp(im * pi/4.0);//std::exp(im*pi/4.0);// exp(2*pi*array[i])+(double)n;
  //   // printf("%f,", real(nums_complex[i]));
  //   // printf("%f,", imag(nums_complex[i]));
  //   // printf("expend\n");
  //   // printf("world_rank:  with number %f\n",nums[i] );
  // }

  // printf("array: %d\n",n );
  // for(int i=0;i<n_qft_point;i++)
  // {
  //   printf("%f,",array[i] );
  // }
  // printf("\n");

  std::vector < ublas::vector<std::complex<double> > > v;
int kkk=2;
  for (int i = n_qft_point-1; i >= 0; i-=1)
  {
    ublas::vector<std::complex<double> > v1(kkk);

    v1(0)=1.0;
    v1(1)=std::exp(2.0*im * pi * array[i]);
    v.push_back(v1);
  }
 // printf("expstart: %d\n",n);
 //  for (int i=0; i<n_qft_point; i++) 
 //    {
 //      printf("%f,%f\n",real(v[i](0) ),imag(v[i](0) ) );
 //      printf("%f,%f\n",real(v[i](1) ),imag(v[i](1) ) );
 //    }

// // for loop for tensor
//     ublas::vector<std::complex<double> > v_current = v.back();
//     v.pop_back();
//   for (int i = 0; i < n_qft_point; ++i)
//   {

//     ublas::vector<std::complex<double> > v_sec_current = v.back();
//     v.pop_back();


//     ublas::matrix<std::complex<double> > m(v_current.size(), v_sec_current.size());
//       m =  outer_prod(v_current, v_sec_current);


//     ublas::vector<std::complex<double> > v_from_m;




//   }






  ublas::matrix<std::complex<double> > m(v[0].size(), v[1].size());
  m =  outer_prod(v[0], v[1]);


  printf("expstart: %d\n",n);
  for(int i = 0;i< m.size1();i++)
  {
    for(int j=0;j<m.size2();j++)
    {
      printf("%f,%f\n",real(m(i,j) ),imag(m(i,j) ) );
    }
  }


  


  // ublas::vector<std::complex<double> >v1(4), v2(2);



  int len_nums = 2*pow(2,n_qft_point);
  double *nums = (double *)malloc(sizeof(double) * len_nums);
  assert(nums != NULL);


  // bringing the nums_complex 1d array to 
  for (int i = 0; i < len_nums_complex; i++) {

    nums[2*i] = real(nums_complex[i]);
    nums[2*i+1]= imag(nums_complex[i]);
    //printf("world_rank: %d , %dth real number %f, img number %f\n",n,i,nums[2*i],nums[2*i+1] );
  }
  return nums;
}
// Computes the average of an array of numbers
double compute_avg(double *array, int num_elements) {
  double sum = 0.f;
  int i;
  for (i = 0; i < num_elements; i++) {
    sum += array[i];
  }
  return sum / num_elements;
}
// Computes the average of an array of numbers
double *compute_neg(double *array, int num_elements,int n) {
  double *nums = (double *)malloc(sizeof(double) * num_elements);
  assert(nums != NULL);
    ublas::vector<std::complex<double> > v1(4), v2(2);
  int i;
  for (i = 0; i < num_elements; i++) {

    nums[i] = array[i]+n;
    //printf("world_rank: %d with number %f\n",n,nums[i] );
  }
  return nums;
}


int main(int argc, char** argv) {






   // std::string binary = std::bitset<8>(pow(2,3)).to_string(); //to binary
   //  std::cout<<binary[0]<<"\n";

   //  unsigned long decimal = std::bitset<8>(binary).to_ulong();
   //  std::cout<<decimal<<"\n";






  std::complex<double> I(0, 1); 
  pi = 4 * atan(1.0);
  // std::complex<double> e(-2*pi/2,0);
  // std::complex<double> z = exp(I*e);
  // double x = 2;
  // std::cout << real(z*2.0) <<'\n';
  ublas::vector<std::complex<double> > v1(4), v2(2);
  // for (unsigned i = 0; i < 4; ++i)
  //     v1 (i) = i;
  //  for (unsigned i = 0; i < 2; ++i)
  //     v2 (i) = i;   

  // ublas::matrix<double> m(v1.size(), v2.size());
  // std::cout 
  //           << v1 << '\n'
  //           << v2 << '\n'
  //           <<  outer_prod(v1, v2)<< '\n';
  int num_point_qft = 2;//atoi(argv[1]);

  // Seed the random number generator to get different results each time


  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  //printf("world_size: %d\n",world_size );

  // Create a random array of elements on the root process. Its total
  // size will be the number of elements per process times the number
  // of processes


        int num_states = pow(2,num_point_qft);
  double *rand_nums = NULL;
  if (world_rank == 0) {



        double final_const = 1/(pow(2,num_point_qft/2));



        // int *decimal_bits = (int *)malloc(sizeof(int) * num_states);
        std::string binary_bits[num_states];
        for (int i = 0; i < num_states; ++i)
        {
          binary_bits[i] = std::bitset<2>(i).to_string(); // thisssss <3> cannot be replaced with num_point_qft
          //std::cout<<binary_bits[i][0]<<"\n";
        }
        double phases[num_states][num_point_qft];

        for (int i = 0; i < num_states; ++i)
        {
          
          for (int j = 0; j < num_point_qft; ++j)
          {
            int start_checking_point = num_point_qft - 1 -j;
            for(int k=0; k < num_point_qft; ++k)
            {
              if(start_checking_point<num_point_qft && binary_bits[i][j]=='1')
              {
                phases[i][start_checking_point]+=pow(2,-(k+1));
              }
              start_checking_point+=1;
            }
          }
        }

        double oneline_phases[num_states*num_point_qft];  
        int k=0;
        for (int i = 0; i < num_states; ++i)
        {
          for(int j=0; j< num_point_qft;++j)
          {
            //printf("%f\t",phases[i][j] );
            oneline_phases[k]=phases[i][j];
            k++;
          }
          //printf("\n");
        }
        // printf("start oneline \n");

        // for (int i = 0; i < num_point_qft*num_states; ++i)
        // {
        //   printf("%f,",oneline_phases[i] );
        // }
        // printf("end oneline \n");


    // rand_nums = create_nums(num_point_qft * world_size,world_rank);
    rand_nums = create_onelinephases(oneline_phases, num_point_qft*num_states);
  }

  // For each process, create a buffer that will hold a subset of the entire
  // array
  double *sub_rand_nums = (double *)malloc(sizeof(double) * num_point_qft);
  assert(sub_rand_nums != NULL);

  // Scatter the random numbers from the root process to all processes in
  // the MPI world
  MPI_Scatter(rand_nums, num_point_qft, MPI_DOUBLE, sub_rand_nums,
              num_point_qft, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Compute the average of your subset
  double *rec_indi_nums = NULL;
  rec_indi_nums = compute_qft(sub_rand_nums, num_point_qft,world_rank);
 
  // Gather all partial averages down to the root process
  double *rec_gather_nums = NULL;
  if (world_rank == 0) {
    rec_gather_nums = (double *)malloc(sizeof(double) * num_states*num_states*2);
    assert(rec_gather_nums != NULL);
  }
  MPI_Gather(rec_indi_nums, num_states*2 , MPI_DOUBLE, rec_gather_nums, num_states*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Now that we have all of the partial averages on the root, compute the
  // total average of all numbers. Since we are assuming each process computed
  // an average across an equal amount of elements, this computation will
  // produce the correct answer.
  if (world_rank == 0) {
    //double avg = compute_avg(sub_avgs, world_size);
    // printf("Avg of all elements is %f\n", avg);
    // // Compute the average across the original data for comparison
    // double original_data_avg =
    //   compute_avg(rand_nums, num_point_qft * world_size);
    // printf("Avg computed across original data is %f\n", original_data_avg);

  // for (int i = 0; i < num_states; i++) 
  // {
  //   printf("%dth state:\n", i);
  //   for (int j = 0; j  < num_states; j++)
  //   {
  //     printf("\t%dth:( ",j );

  //       printf("%f , %f )\n",rec_gather_nums[ (2*num_states*i)+(2*j) ],rec_gather_nums[ (2*num_states*i)+(2*j+1)] );

  //   }
  // }

  // ublas::vector<std::complex<double> > v1(4), v2(2);
  //  for (double i = 0; i < 4; ++i)
  //     v1 (i) = exp(I*i);
  //  for (double i = 0; i < 2; ++i)
  //     v2 (i) = i;   

  // ublas::matrix<double> m(v1.size(), v2.size());
  // std::cout 
  //           << v1 << '\n'
  //           << v2 << '\n'
  //           <<  outer_prod(v1, v2)<< '\n';
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
