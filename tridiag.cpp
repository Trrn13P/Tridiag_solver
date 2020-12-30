#include "tridiag.hpp"

#include <armadillo>
using namespace arma;


//General forward solver
void tridiag::forward_solver(vec &v_vec, vec &g_vec){
  //Setting inital values for the tilde vectors
  b_tilde(0) = b(0); c_tilde(0) = c(0); g_tilde(0) = g_vec(0);


  //precalc
  float d;
  //Solving the forward algorythm
  for(int i=1;i<n;i++){
    d = b_tilde(i-1)*1./a(i-1);

    b_tilde(i) = b(i)*d - c_tilde(i-1);
    c_tilde(i) = c(i)*d;
    g_tilde(i) = g_vec(i)*d - g_tilde(i-1);

  }

}


//General backward_solver
void tridiag::backward_solver(vec &v_vec, vec &g_vec){
  //Setting up the endpoint
  v_vec(n) = g_tilde(n-1)*1./b_tilde(n-1);


  //Solving the backward soloution
  for(int i=2;i<n+1;i++){
    int j = (n-i);
    v_vec(j+1) = (g_tilde(j) - c_tilde(j) *v_vec(j+2))*1./b_tilde(j);
  }

}

void tridiag::solve(vec &v_vec, vec &g_vec){
  forward_solver(v_vec, g_vec);
  backward_solver(v_vec, g_vec);
}
