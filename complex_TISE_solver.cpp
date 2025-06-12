#include <Eigen> // Eigen for now, later maybe LAPACK
#include <iostream>
const double h_bar = 1.0;
const double m = 1.0;

std::complex<double> V(double x)
{
	return std::complex<double>(0.0, x*x*x);
}


int main()
{
	const double x_min = -5;
	const double x_max = 5;
	
	const int N = 5;
	const int N_evals = N; //edge cases
	
	const double dx = (x_max - x_min)/(N-1);
	float x[N_evals]; x[0] = x_min;
	for (int i=1; i<N_evals; ++i)
	{
		x[i] = x[i-1]+dx;
	}

	Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(N_evals, N_evals); //hamiltonian

	const std::complex<double> KE = - (h_bar*h_bar)/(2*m*dx*dx); //finite difference KE with the (dx*dx) from the derivative approximation
	//fill hamiltonian.
	H(0, 0) = V(x_min)-2.0*KE; H(0, 1) = KE;
	H(N_evals-1, N_evals-1) = V(x_max-2*dx) - 2.0*KE; H(N_evals-1, N_evals-2) = KE;
	for (int i=1; i<N_evals-1; ++i)
	{
		H(i, i) = V(x[i])-2.0*KE;
		H(i, i-1) = KE; H(i, i+1) = KE;
	}
	std::cout << H << std::endl;
	std::cout << "----------------" << std::endl;
	
	//solve eigval
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es; // solve complex self-adjoint matrix eigenvector problem
	es.compute(H);
	//for (int i=0; i<N; ++i)

	std::cout << es.eigenvectors() << std::endl;

	/*
	std::cout << "START" << std::endl;
	const int maxrets = 5;
	int ret_count = 0;
	for (int j=0; j<N; ++j)
	{
		if (std::abs(es.eigenvalues()(j).imag()) < 0.1 && ret_count<maxrets)
		{
			for (int i=0; i<N; ++i)
			{
				std::cout << "x = " << x[i] << ", psi = " << es.eigenvectors()(i, j) << std::endl;
			}
			std::cout << "-" << std::endl;
			++ret_count;
		}
	}
	std::cout << "STOP" << std::endl;
	*/
	
	std::cout << "-----------------" << std::endl;
	std::cout << es.eigenvalues() << std::endl;
	return 0;
}	
		
