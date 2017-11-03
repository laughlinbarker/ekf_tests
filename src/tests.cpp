/*
This is an implementation of real-time statistical tests found in Chapter 5, Bar-Shalom, 2001, for examining EKF consistancy.
See page 237 for details.

Laughlin Barker, Nov. 2017, Dynamical Systems and Control Laboratory, JHU
laughlinbarker@gmail.com
*/

#ifndef EKF_TESTS_
#define EKF_TESTS_

//Eigen for Linear Algebra
#include <Eigen/Dense>

//boost for chi-sq distribution
#include <boost/math/distributions/chi_squared.hpp>

#include <iostream>


class EKFConsistancyTests
{
public:

	EKFConsistancyTests(int N, double alpha)
	{
    N = N;
	}

  //Compute single innovation squared statistic
  double computeInnovSqStat(Eigen::VectorXd &v, Eigen::MatrixXd &S)
  {
    //ensure v and S are appropriate size
    if ((S.rows() == S.cols()) && (S.rows() == v.size()))
      return v.transpose() * S.inverse() * v;
    else
      return -1;
  }

  bool passChiSqHypothesis(int N, double alpha, double val)
  {
    //initialize chi-sq distribution
    chi_squared dist(N);

    //compute upper and lower quintile bounds (using half of distribution)
    double upperQ = quantile(complement(dist, alpha / 2));
    double lowerQ = quantile(dist, alpha / 2);

    using boost::math::chi_squared;
    using boost::math::quantile;
    using boost::math::complement;

  }
private:

	//number of samples over which to evaluate the statistic
	int N;
	//tail probability 
	double alpha;

};

int main()
{

  EKFConsistancyTests test(2,0.05);

  Eigen::VectorXd v(3);
  v << 1, 1, 1;
  Eigen::MatrixXd S(3,3);
  S = Eigen::MatrixXd::Identity(3,3);

  std::cout <<"S: " << S<< std::endl;

  std::cout<< "Quadratic: "<< test.computeInnovSqStat(v,S) << std::endl;
}

#endif