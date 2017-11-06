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

//for chi-sq distribution, and rolling mean
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>

#include <iostream>


using namespace boost::accumulators;
class EKFConsistancyTests
{

public:

	EKFConsistancyTests(int N, double alpha, int nz):
  acc(tag::rolling_window::window_size = N)
	{
    using boost::math::chi_squared;
    using boost::math::quantile;
    using boost::math::complement;

    //store input vars
    N = N;
    alpha = alpha;
    nz = nz;

    //compute appropriate interval for Chi-Squared test
    chi_squared dist(N * nz);
    upperQ = (quantile(complement(dist, alpha / 2)) )/ N;
    lowerQ = (quantile(dist, alpha / 2) ) / N;

    std::cout << "Interval: ["<< lowerQ << "," << upperQ << "] \n";

    //set test counters to zero
    iSqTest, iAutoCorr = 0;
  }

  //Compute single innovation squared statistic
  double computeInnovSqStat(Eigen::VectorXd &v, Eigen::MatrixXd &S)
  {
    //ensure v and S are of appropriate size
    if ((S.rows() == nz) && (S.cols() == nz) && (v.size() == nz))
      return v.transpose() * S.inverse() * v;
    else
    {
      std::cout << "ERROR: input vec. or matrix size don't agree with size given in initialization..." << std::endl;
      return -1;
    }
  }

  bool chiSqHypothesisTest(double value)
  {
    if ( (value < lowerQ) || (value > upperQ) )
    {
      //value outisde of probability regions for acceptance
      return false;
    }
    else
    {
      //value inside interval for acceptance
      return true;
    }
  }

  int checkInnovationStatistic(Eigen::VectorXd &v, Eigen::MatrixXd &S)
  {
    //filter not ready, update count and 
    if (iSqTest < N)
    {
      iSqTest++;
    }

    double innovSqStat = computeInnovSqStat(v,S);

    if (innovSqStat >= 0)
    { 
      //add most recent innovation squared statistic value to the rolling-average accumulator
       acc(innovSqStat);
   
       double mean = rolling_mean(acc);
   
       //check to see if filter is fully initialized
       if (iSqTest >= N)
       {
         innovationSqTestReady = true;
         if ((mean > lowerQ) && (mean < upperQ))
         {
           //test passed, filter appears 
           return 1;
         }
         else
         {
           return -1;
         }
       }
   
       //filter not ready yet, return 0
       else if
       {
         return 0;
       }
    }
    else
      //error, wrong size vec or matrix
      return -1;
  }

private:

	//number of samples over which to evaluate the statistic
	int N;

  //number of times we've called the two filters
  int iSqTest, iAutoCorr;
	//tail probability 
	double alpha;
  //size of innovation vec.
  int nz;

  //Chi-sq bounds
  double upperQ, lowerQ;

  //flag as to wether the filter has been called sufficient (N) number of times
  //before we can trust the test's output
  bool innovationSqTestReady;
  bool timeAvgAutoCorrelationTestReady;

  //rolling average accumulator
  accumulator_set<double, stats<tag::rolling_mean> > acc;

};

int main()
{

  //number of samples over which we wish to compute the time-
  //average normalized innovation-squared statistic
  int N = 50;
  // Confidance interval for statistica test (1-alpha)
  double alpha = 0.05;
  // Size of innovation vector
  int z = 2;

  EKFConsistancyTests test(N,alpha,z);

  Eigen::VectorXd v(3);
  v << 1, 1, 1;
  Eigen::MatrixXd S(3,3);
  S = Eigen::MatrixXd::Identity(3,3);

}


#endif