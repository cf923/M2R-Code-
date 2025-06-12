#include <Eigen> // Eigen for now, later maybe LAPACK
#include <iostream>
//note: precision overestimation - way too much stuff is e.g. double when it could be float, but this can be fixed later

// chuck constants in the initialized data segment for now, could be preprocessor directives (should be)
// if you set these to be physically meaningful (e.g. correct values of h_bar), plugging in x on amstrong lengths and electron volt potentials will give dimensionally correct outputs.
const double h_bar = 1.0;
const double m = 1.0;

double V(double x){
	/*
	// well, well, well...
	if (x>=-2.0 && x<= 2.0){
		return 0.0;
	} else {
		return 10000.0;
	}
	*/
	//simple harmonic oscillator:
	//const double k = 1.0;
	
	//return x*x;
        if (x >= -1 && x <= 1)
        {
                return -3.0;
        }
        else if (x >= -3 &&  x < -1)
        {
                return -10.0;
        }
        else if (x > 1 and x <= 3)
        {
                return -5.0;
        }
        else
        {
                return 10000.0;
        }

}
int main(){
	const double x_min = -5;
	const double x_max = 5;
	
	const int N = 250;
	const int N_evals = N-2; //stupid edge cases
	
	const double dx = (x_max - x_min)/(N-1);
	float x[N_evals]; x[0] = x_min;
	for (int i=1; i<N_evals; ++i){
		x[i] = x[i-1]+dx;
	}

	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N_evals, N_evals); //hamiltonian

	const double KE = - (h_bar*h_bar)/(2*m*dx*dx); //finite difference KE with the (dx*dx) from the derivative approximation
	//fill hamiltonian. It looks awful but it's to remove conditionals while dealing with edge cases. Faster.
	H(0, 0) = V(x_min)-2.0*KE; H(0, 1) = KE;
	H(N_evals-1, N_evals-1) = V(x_max-2*dx) - 2.0*KE; H(N_evals-1, N_evals-2) = KE;
	for (int i=1; i<N_evals-1; ++i){
		H(i, i) = V(x[i])-2.0*KE;
		H(i, i-1) = KE; H(i, i+1) = KE;
	}
        for (int i=0; i<N_evals; ++i){
                for (int j =0; j<N_evals; ++j){
                        std::cout << " " <<  H(i, j);
                }
                std::cout << "\n";
	}
	
	//solve eigval
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H); // i am not writing an eigenvalue solver myself
	// this is really bad but it works. the gnuplot implementation is much nicer (i'm lying)
  	for (int j=0; j<10; ++j){ // the upper limit in this gives the number of states returned. it gets hard to see with too high a number, i suggest below 10.
		//again, change these with some boundary conditions
   		std::cout << "x = " << x_min << ", psi = 0.0" << std::endl;
		for (int i = 0; i < N_evals; ++i) {
        		// Eigenvector elements correspond to psi_1, psi_2, ..., psi_{N-2}
        		std::cout << "x = " << x[i] << ", psi = " << es.eigenvectors()(i, j) << std::endl;
    		}
    		std::cout << "x = " << x_max << ", psi = 0.0" << std::endl;
		std::cout << "-" << std::endl; // indicates new vector to the python script that parses this
	}
	return 0;
}	
		
