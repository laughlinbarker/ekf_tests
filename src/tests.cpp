/*
This is an implementation of real-time statistical tests found in Chapter 5, Bar-Shalom, 2001, for examining EKF consistancy.
See page 237 for details.

Laughlin Barker, Nov. 2017, Dynamical Systems and Control Laboratory, JHU
laughlinbarker@gmail.com
*/

#ifndef EKF_TESTS_
#define EKF_TESTS_

// //Eigen for Linear Algebra
// #include <Eigen/Core>

// //boost for chi-sq distribution
// #include <boost/math/distributions/chi_squared.hpp>

// //misc
// #include <iostream>


// class EKFConsistancyTests
// {
// public:

// 	EKFConsistancyTests(int N, float alpha)
// 	{

// 	}


// private:
// 	//measurement innovation vector (v = z - z_predicted)
// 	Eigen::VectorXf v;
// 	//measurement prediction covariance matrix (S = HPH' + R) 
// 	Eigen::MatrixXf S;

// 	//number of samples over which to evaluate the statistic
// 	int N;
// 	//tail probability 
// 	double alpha;

// 	boost::math::chi_squared dist;

// 	//Compute single innovation squared statistic
// 	double computeInnovSqStat(Eigen::VectorXf v, Eigen::MatrixXf S)
// 	{
// 		//ensure v and S are appropriate size
// 		if (S.rows() == S.cols()) && (S.rows() == v.size())
// 			return v.transpose() * S.inverse() * v;
// 		else
// 			return -1;
// 	}


// };
// Copyright John Maddock 2006, 2007
// Copyright Paul A. Bristow 2010

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
using std::cout; using std::endl;
using std::left; using std::fixed; using std::right; using std::scientific;
#include <iomanip>
using std::setw;
using std::setprecision;

#include <boost/math/distributions/chi_squared.hpp>

void confidence_limits_on_std_deviation(
        double Sd,    // Sample Standard Deviation
        unsigned N)   // Sample size
{
   // Calculate confidence intervals for the standard deviation.
   // For example if we set the confidence limit to
   // 0.95, we know that if we repeat the sampling
   // 100 times, then we expect that the true standard deviation
   // will be between out limits on 95 occations.
   // Note: this is not the same as saying a 95%
   // confidence interval means that there is a 95%
   // probability that the interval contains the true standard deviation.
   // The interval computed from a given sample either
   // contains the true standard deviation or it does not.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm

   // using namespace boost::math; // potential name ambiguity with std <random>
   using boost::math::chi_squared;
   using boost::math::quantile;
   using boost::math::complement;

   // Print out general info:
   cout <<
      "________________________________________________\n"
      "2-Sided Confidence Limits For Standard Deviation\n"
      "________________________________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Number of Observations" << "=  " << N << "\n";
   cout << setw(40) << left << "Standard Deviation" << "=  " << Sd << "\n";
   //
   // Define a table of significance/risk levels:
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Start by declaring the distribution we'll need:
   chi_squared dist(N - 1);
   //
   // Print table header:
   //
   cout << "\n\n"
           "_____________________________________________\n"
           "Confidence          Lower          Upper\n"
           " Value (%)          Limit          Limit\n"
           "_____________________________________________\n";
   //
   // Now print out the data for the table rows.
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // Calculate limits:
      double lower_limit = sqrt((N - 1) * Sd * Sd / quantile(complement(dist, alpha[i] / 2)));
      double upper_limit = sqrt((N - 1) * Sd * Sd / quantile(dist, alpha[i] / 2));
      // Print Limits:
      cout << fixed << setprecision(5) << setw(15) << right << lower_limit;
      cout << fixed << setprecision(5) << setw(15) << right << upper_limit << endl;
   }
   cout << endl;
} // void confidence_limits_on_std_deviation

void confidence_limits_on_std_deviation_alpha(
        double Sd,    // Sample Standard Deviation
        double alpha  // confidence
        )
{  // Calculate confidence intervals for the standard deviation.
   // for the alpha parameter, for a range number of observations,
   // from a mere 2 up to a million.
   // O. L. Davies, Statistical Methods in Research and Production, ISBN 0 05 002437 X,
   // 4.33 Page 68, Table H, pp 452 459.

   //   using namespace std;
   // using namespace boost::math;
   using boost::math::chi_squared;
   using boost::math::quantile;
   using boost::math::complement;

   // Define a table of numbers of observations:
   unsigned int obs[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40 , 50, 60, 100, 120, 1000, 10000, 50000, 100000, 1000000};

   cout <<   // Print out heading:
      "________________________________________________\n"
      "2-Sided Confidence Limits For Standard Deviation\n"
      "________________________________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Confidence level (two-sided) " << "=  " << alpha << "\n";
   cout << setw(40) << left << "Standard Deviation" << "=  " << Sd << "\n";

   cout << "\n\n"      // Print table header:
            "_____________________________________________\n"
           "Observations        Lower          Upper\n"
           "                    Limit          Limit\n"
           "_____________________________________________\n";
    for(unsigned i = 0; i < sizeof(obs)/sizeof(obs[0]); ++i)
   {
     unsigned int N = obs[i]; // Observations
     // Start by declaring the distribution with the appropriate :
     chi_squared dist(N - 1);

     // Now print out the data for the table row.
      cout << fixed << setprecision(3) << setw(10) << right << N;
      // Calculate limits: (alpha /2 because it is a two-sided (upper and lower limit) test.
      double lower_limit = sqrt((N - 1) * Sd * Sd / quantile(complement(dist, alpha / 2)));
      double upper_limit = sqrt((N - 1) * Sd * Sd / quantile(dist, alpha / 2));
      // Print Limits:
      cout << fixed << setprecision(4) << setw(15) << right << lower_limit;
      cout << fixed << setprecision(4) << setw(15) << right << upper_limit << endl;
   }
   cout << endl;
}// void confidence_limits_on_std_deviation_alpha

void chi_squared_test(
       double Sd,     // Sample std deviation
       double D,      // True std deviation
       unsigned N,    // Sample size
       double alpha)  // Significance level
{
   //
   // A Chi Squared test applied to a single set of data.
   // We are testing the null hypothesis that the true
   // standard deviation of the sample is D, and that any variation is down
   // to chance.  We can also test the alternative hypothesis
   // that any difference is not down to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm
   //
   // using namespace boost::math;
   using boost::math::chi_squared;
   using boost::math::quantile;
   using boost::math::complement;
   using boost::math::cdf;

   // Print header:
   cout <<
      "______________________________________________\n"
      "Chi Squared test for sample standard deviation\n"
      "______________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(55) << left << "Number of Observations" << "=  " << N << "\n";
   cout << setw(55) << left << "Sample Standard Deviation" << "=  " << Sd << "\n";
   cout << setw(55) << left << "Expected True Standard Deviation" << "=  " << D << "\n\n";
   //
   // Now we can calculate and output some stats:
   //
   // test-statistic:
   double t_stat = (N - 1) * (Sd / D) * (Sd / D);
   cout << setw(55) << left << "Test Statistic" << "=  " << t_stat << "\n";
   //
   // Finally define our distribution, and get the probability:
   //
   chi_squared dist(N - 1);
   double p = cdf(dist, t_stat);
   cout << setw(55) << left << "CDF of test statistic: " << "=  "
      << setprecision(3) << scientific << p << "\n";
   double ucv = quantile(complement(dist, alpha));
   double ucv2 = quantile(complement(dist, alpha / 2));
   double lcv = quantile(dist, alpha);
   double lcv2 = quantile(dist, alpha / 2);
   cout << setw(55) << left << "Upper Critical Value at alpha: " << "=  "
      << setprecision(3) << scientific << ucv << "\n";
   cout << setw(55) << left << "Upper Critical Value at alpha/2: " << "=  "
      << setprecision(3) << scientific << ucv2 << "\n";
   cout << setw(55) << left << "Lower Critical Value at alpha: " << "=  "
      << setprecision(3) << scientific << lcv << "\n";
   cout << setw(55) << left << "Lower Critical Value at alpha/2: " << "=  "
      << setprecision(3) << scientific << lcv2 << "\n\n";
   //
   // Finally print out results of alternative hypothesis:
   //
   cout << setw(55) << left <<
      "Results for Alternative Hypothesis and alpha" << "=  "
      << setprecision(4) << fixed << alpha << "\n\n";
   cout << "Alternative Hypothesis              Conclusion\n";
   cout << "Standard Deviation != " << setprecision(3) << fixed << D << "            ";
   if((ucv2 < t_stat) || (lcv2 > t_stat))
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << "Standard Deviation  < " << setprecision(3) << fixed << D << "            ";
   if(lcv > t_stat)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << "Standard Deviation  > " << setprecision(3) << fixed << D << "            ";
   if(ucv < t_stat)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
} // void chi_squared_test

void chi_squared_sample_sized(
        double diff,      // difference from variance to detect
        double variance)  // true variance
{
   using namespace std;
   // using boost::math;
   using boost::math::chi_squared;
   using boost::math::quantile;
   using boost::math::complement;
   using boost::math::cdf;

   try
   {
   cout <<   // Print out general info:
     "_____________________________________________________________\n"
      "Estimated sample sizes required for various confidence levels\n"
      "_____________________________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(40) << left << "True Variance" << "=  " << variance << "\n";
   cout << setw(40) << left << "Difference to detect" << "=  " << diff << "\n";
   //
   // Define a table of significance levels:
   //
   double alpha[] = { 0.5, 0.33333333333333333333333, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Print table header:
   //
   cout << "\n\n"
           "_______________________________________________________________\n"
           "Confidence       Estimated          Estimated\n"
           " Value (%)      Sample Size        Sample Size\n"
           "                (lower one-         (upper one-\n"
           "                 sided test)        sided test)\n"
           "_______________________________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // Calculate df for a lower single-sided test:
      double df = chi_squared::find_degrees_of_freedom(
         -diff, alpha[i], alpha[i], variance);
      // Convert to integral sample size (df is a floating point value in this implementation):
      double size = ceil(df) + 1;
      // Print size:
      cout << fixed << setprecision(0) << setw(16) << right << size;
      // Calculate df for an upper single-sided test:
      df = chi_squared::find_degrees_of_freedom(
         diff, alpha[i], alpha[i], variance);
      // Convert to integral sample size:
      size = ceil(df) + 1;
      // Print size:
      cout << fixed << setprecision(0) << setw(16) << right << size << endl;
   }
   cout << endl;
   }
  catch(const std::exception& e)
  { // Always useful to include try & catch blocks because default policies
    // are to throw exceptions on arguments that cause errors like underflow, overflow.
    // Lacking try & catch blocks, the program will abort without a message below,
    // which may give some helpful clues as to the cause of the exception.
    std::cout <<
      "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
  }
} // chi_squared_sample_sized

int main()
{
   // Run tests for Gear data
   // see http://www.itl.nist.gov/div898/handbook/eda/section3/eda3581.htm
   // Tests measurements of gear diameter.
   //
   confidence_limits_on_std_deviation(0.6278908E-02, 100);
   chi_squared_test(0.6278908E-02, 0.1, 100, 0.05);
   chi_squared_sample_sized(0.1 - 0.6278908E-02, 0.1);
   //
   // Run tests for silicon wafer fabrication data.
   // see http://www.itl.nist.gov/div898/handbook/prc/section2/prc23.htm
   // A supplier of 100 ohm.cm silicon wafers claims that his fabrication
   // process can produce wafers with sufficient consistency so that the
   // standard deviation of resistivity for the lot does not exceed
   // 10 ohm.cm. A sample of N = 10 wafers taken from the lot has a
   // standard deviation of 13.97 ohm.cm
   //
   confidence_limits_on_std_deviation(13.97, 10);
   chi_squared_test(13.97, 10.0, 10, 0.05);
   chi_squared_sample_sized(13.97 * 13.97 - 100, 100);
   chi_squared_sample_sized(55, 100);
   chi_squared_sample_sized(1, 100);

   // List confidence interval multipliers for standard deviation
   // for a range of numbers of observations from 2 to a million,
   // and for a few alpha values, 0.1, 0.05, 0.01 for condfidences 90, 95, 99 %
   confidence_limits_on_std_deviation_alpha(1., 0.1);
   confidence_limits_on_std_deviation_alpha(1., 0.05);
   confidence_limits_on_std_deviation_alpha(1., 0.01);

   return 0;
}
#endif