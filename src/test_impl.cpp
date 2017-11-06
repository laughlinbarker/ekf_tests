#include <ekf_tests/tests.cpp>

#include <random>

int main()
{

  //number of samples over which we wish to compute the time-
  //average normalized innovation-squared statistic
  int N = 20;
  // Confidance interval for statistica test (1-alpha)
  double alpha = 0.05;
  // Size of innovation vector
  int z = 2;

  EKFConsistancyTests test(N,alpha,z);

  //create fake covariance matrix, which corresponds to 3 IID RV with variance = 1
  Eigen::VectorXd v(z);
  Eigen::MatrixXd S(z,z);
  S = Eigen::MatrixXd::Identity(z,z);

  //setup RNG with N(0,1) distribution
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0,1);

  int result = 0;

  //iterate through test shoving new vectors into filter each time...
  for (int i=0; i<200; i++){

  	for (int ii = 0; ii<v.size(); ii++){
  		v(ii) = distribution(generator);
  	}

  	//add some bias half way through
  	if (i > 100)
  	{
  		v(0) += 1;
  		v(1) += 1;
  	}

  	result = test.checkInnovationStatistic(v,S);

  	std::cout << "Iteration: " << i << std::endl;

  	if (result == 0)
  		std::cout << "Filter not initialized" << std::endl;
  	else if (result  == 1)
	  	std::cout <<  "Appears consistant!" << std::endl;
  	else
  		std::cout << "Appears to have lost constancy..." << std::endl;
  }

}