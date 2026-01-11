#include "tools.hpp"
#include <mpi.h>

int main(int argc, char**argv) {
  MPI_Init(&argc, &argv); 
  int rank, world_size; 

  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &world_size); 

  if(rank == 0){
   std::cout<<"World size:"<<world_size<<std::endl;
  }

  // Things to the matriz
  std:: vector < double > conini = {0, 100};
  size_t dim{9}; // by default

  if (argc > 1) {
    dim = std:: atoi(argv[1]); // for leonardo 12000
  }

  size_t ite{1000}; // leonardo 1000
  size_t printInterval{200};

  // Set the things to each rank
  size_t dim_local = dim/world_size; 
  size_t rest = dim % world_size; 
  size_t start_row; 

  if (rank < rest) {
    dim_local += 1; 
    start_row = rank*dim_local; 
  } else {
    start_row = rank*dim_local+rest; 
  }

  // Create the mesh class element
  CMesh < double > Matrix(dim, dim_local, start_row, conini, rank, world_size); 
  // Matrix.print_in_serial(); 
  //Matrix.print_in_parallel(); 
  //

  CSolver < double > solver(Matrix); 
  {
    Parallel_Timer t("Total_Time"); 
    solver.jacobi(Matrix, ite, printInterval); 
  }
  if (rank == 0) {
    std:: cout << "Remember comp_time = comp_time-comm_time "<<std:: endl; 
  }
  //Matrix.print_in_parallel(); // Matrix after all the iterations
  //
  //  print the times
  Parallel_Timer:: gather_timing_data(MPI_COMM_WORLD, 0); 

  MPI_Finalize(); 
  return 0; 
}
