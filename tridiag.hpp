#include <armadillo>
using namespace arma;

class tridiag {

  private:
    //Setting up integers
    int n;
    //Setting up vectors
    vec a, b, c, g_vec;
    vec b_tilde, c_tilde, g_tilde;
    vec v_vec, x_vec, u_vec;





    void Initialize(int n_, vec a_, vec b_, vec c_){
      a = a_; b = b_; c = c_;

      v_vec = zeros(n+1);


      //setting up the tilde vectors
      b_tilde = zeros(n+1);
      c_tilde = zeros(n+1);
      g_tilde = zeros(n+1);


    }



public:
  void forward_solver(vec &v_vec, vec &g_vec);
  void backward_solver(vec &v_vec, vec &g_vec);
  void solve(vec &v_vec, vec &g_vec);



  tridiag(int n_, vec a_, vec b_, vec c_){
    tridiag(n_, a_, b_, c_);
  }

};
