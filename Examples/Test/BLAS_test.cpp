#include <iostream>
#include <vector>

using namespace std;

// C++ calls functions in a different way, so you need to change specify that this is a C/FORTRAN function call

extern "C"{
// FORTRAN adds _ after all the function names
// and all variables are called by reference
double ddot_( const int *N, const double *a, const int *inca, const double *b, const int *incb );
}

double ddot( int N, double *a, int inca, double *b, int incb ){
  return ddot_( &N, a, &inca, b, &incb );
};

double ddot_vec( std::vector<double> v, std::vector<double> w){
  int N = v.size();
  int i( 1 );
  return ddot_( &N, &*v.begin(), &i, &w[0], &i);
}

int main( int argc, char** argv ){
  // you can define the arrays in one of two ways
  // on the heap
  double *a = new double[3];
  a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
  // on the stack
  double b[3] = { 4.0, 5.0, 6.0 };

  double dot_product = ddot( 3, a, 1, b, 1 );
  cout <<"The dot product is: " <<  dot_product << endl;

  vector<double> v(3), w(3);
  v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;
  w[0] = 4.0; w[1] = 5.0; w[2] = 6.0;

  double dot_prod = ddot_vec( v, w);
  cout << "v.w = " << dot_prod << endl;

  return 0;
};
