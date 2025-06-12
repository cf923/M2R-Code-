#include <Eigen> // Eigen for now, later maybe LAPACK
#include <iostream>
#include <typeinfo>

const double h_bar = 1.0;
const double m = 1.0;

const std::complex<double> imagi(0.0, 1.0);

std::complex<double> func(std::complex<double> x)
{
	return imagi*x*x*x;
	
	/*
	ENERGY WELL - THIS WORKS
	if (x>=-2 && x<=2)
	{
		return 0;
	} else
	{
		return 10000;
	}
	*/
	
	//return std::complex<double>(-3.0, 2.0);
}


int main()
{
	const double x_min = -5;
	const double x_max = 5;
	
	const int N = 400;
	const int N_evals = N; //edge cases
	const double delta_x = (x_max - x_min)/(N-1);
	std::cout << delta_x << std::endl;
	const std::complex<double> dx(delta_x, 0.0);
	std::cout << dx << std::endl;
	std::complex<double> x[N_evals];
	x[0].real(x_min);
	x[0].imag(0.0);
	//x[0] = (x_min, 0);
	for (int i=1; i<N_evals; ++i)
	{
		std::cout << "x at " << i << ": " << std::endl;
		x[i] = x[i-1]+dx;
	}
	for (int i=0; i<N; ++i)
	{
		std::cout << x[i] << std::endl;
	}

	Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(N_evals, N_evals); //hamiltonian

	const std::complex<double> KE(-(h_bar*h_bar)/(2*m*dx.real()*dx.real()), 0.0); //finite difference KE with the (dx*dx) from the derivative approximation
	std::cout << "KE = " << KE << std::endl;
	//fill hamiltonian.
	//H(0, 0) = V(x_min)-2.0*KE; H(0, 1) = KE;
	//H(N_evals-1, N_evals-1) = V(x_max-2*dx) - 2.0*KE; H(N_evals-1, N_evals-2) = KE;
	
	H(0, 0) = func(x_min)+2.0*KE;
	H(0, 1) = -1.0*KE;
	//H(0, 2) = 4.0*KE;
	//H(0, 3) = -KE;

	H(N_evals-1, N_evals-1) = func(x_max)+2.0*KE;
	H(N_evals-1, N_evals-2) = -1.0*KE;
	//H(N_evals-1, N_evals-3) = 4.0*KE;
	//H(N_evals-1, N_evals-4) = -KE;
	
	//H(1, 1) = V(x[1])-2.0*KE;
	//H(1, 2) = KE;

	//H(N_evals-2, N_evals-2) = V(x[N_evals-2]) - 2.0*KE;
	//H(N_evals-2, N_evals-3) = KE;

	for (int i=1; i<N_evals-1; ++i)
	{
		H(i, i) = func(x[i])-2.0*KE;
		H(i, i-1) = KE;
		H(i, i+1) = KE;
	}
	std::cout << "created H matrix: " << std::endl;
	std::cout << H << std::endl;
	std::cout << "----------------" << std::endl;

	//solve eigval
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces; // solve complex self-adjoint matrix eigenvector problem
	ces.compute(H);

	Eigen::VectorXcd eigenvalues_vec = ces.eigenvalues();
	Eigen::MatrixXcd eigenvectors_mat = ces.eigenvectors();


	Eigen::MatrixXcd lambda_mat = Eigen::MatrixXcd::Zero(H.rows(), H.cols()); // should be replaceable with just NxN
	for (int i=0; i < H.rows(); ++i)
	{
		lambda_mat(i, i) = eigenvalues_vec(i);
	}
	// lambda_mat.diagonal() = eigenvalues_vec; shouuld also work
	
	Eigen::MatrixXcd eigenvectors_inv = eigenvectors_mat.inverse(); // NOTE: may be unstable, especially for large or nearly singular matrices, there are many different methods, if this doesn't work, use another

	Eigen::MatrixXcd H_reconstructed = eigenvectors_mat * lambda_mat  * eigenvectors_inv;

	std::cout << "reconstructed H matrix: " << std::endl;
	std::cout << H_reconstructed << std::endl;

	double error_norm = (H - H_reconstructed).norm();
	std::cout << "Frobenius norm of H - reconstructed H: " << error_norm << std::endl;
	std::cout << "----------------" << std::endl;



	//for (int i=0; i<N; ++i)

	//std::cout << es.eigenvectors() << std::endl;


	std::cout << "START" << std::endl;
	const int maxrets = 15;
	int ret_count = 0;
	for (int j=0; j<N; ++j)
	{
		if (ret_count<maxrets && std::abs(ces.eigenvalues()(j).imag()) < 0.01)
		{
			std::cout << "eigval: " << ces.eigenvalues()(j) << std::endl;
			for (int i=0; i<N; ++i)
			{
				std::cout << "x = " << x[i] << ", psi = " << ces.eigenvectors()(i, j) << std::endl;
			}
			std::cout << "-" << std::endl;
			++ret_count;
		}
	}
	std::cout << "STOP" << std::endl;


	std::cout << "-----------------" << std::endl;
	std::cout << ces.eigenvalues() << std::endl;
	return 0;
}	
