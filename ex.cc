//
// ex.cc - a simple example for exercising the AdaptiveSampler
//         class for generic importance sampling of an unknown
//         population. This example uses a single 5-dimensional
//         Gaussian ellipsoid with correlations between all 5
//         coordinates.
//
// author: richard.t.jones at uconn.edu
// version: march 1, 2017

#include <math.h>
#include <iostream>
#include "AdaptiveSampler.hh"
#include <TRandom.h>

// this example uses the TRandom generator from the CERN ROOT package
TRandom random_generator(0);

void unif01(int n, double *u) { random_generator.RndmArray(n,u); }

int main(int argc, char *argv[])
{
   // set up the adaptive sampler

   const int D=5;
   long int nMC = 1e10;
   AdaptiveSampler sampler(D, &unif01);
   if (argc > 1)
      sampler.restoreState(argv[1]);
   sampler.setVerbosity(2);

   // define the function to be integrated (5D)

   double axis[5][D] = {{ 1, 1, 1, 1, 1},
                        { 1, 1, 1, 1,-4},
                        { 1, 1, 1,-3, 0},
                        { 1, 1,-2, 0, 0},
                        { 1,-1, 0, 0, 0}};
   double mean[D] = {2.0834496,-0.3645377, 1.17464764,-0.398642, 0.6886666};
   double sigma[D] = {0.081, 0.056, 0.021, 0.032, 0.0789};
   double jacobian = 120;
   double gnorm = 1 / jacobian;
   for (int i=0; i < D; ++i)
      gnorm *= sqrt(2*M_PI) * sigma[i];
   double efficiency = gnorm * pow(2, D/2.);

   std::cout << "Problem defined with unbiased efficiency = "
             << efficiency << std::endl
             << "Generating initial sample for adaptation..." 
             << std::endl;

   // perform the Monte Carlo iteration, collect statistics

   int saved_sampler = 0;
   double old_efficiency = efficiency;
   double sum[3] = {0,0,0};
   for (long int i=0; i < nMC; i++) {
      double u[D];
      double w = sampler.sample(u);
      double v[D] = {};
      double gfactor = 0;
      for (int i=0; i < D; ++i) {
         for (int j=0; j < D; ++j) {
            v[i] += axis[i][j] * u[j];
         }
         gfactor += pow((v[i] - mean[i]) / sigma[i], 2);
      }
      double I = exp(-0.5 * gfactor) / gnorm;
      double wI = w * I;
      sum[0] += 1;
      sum[1] += wI;
      sum[2] += wI * wI;
      sampler.feedback(u,wI); // if you want AdaptiveSampler to adapt

      // report progress (or not)

      if (i > 0 && i % 100000000 == 0) {
         double error;
         double result = sampler.getResult(&error);
         std::cout << "result = " << result << " +/- " << error
                   << ", effective sample size " 
                   << sampler.getEfficiency() * sampler.getNsample()
                   << " after " << sampler.getNsample() << " iterations"
                   << std::endl;
         if (sampler.getEfficiency() < efficiency * 0.7 ||
             sampler.getEfficiency() > efficiency * 1.4)
         {
            std::cout << "efficiency expected=" << efficiency 
                      << ", measured=" << sampler.getEfficiency();
            std::cout << ", expand statistics (e)"
                      << " or reset statistics (r)"
                      << " or revert the last update (b)"
                      << " and try again, or accept and go on (a)? ";
            std::string resp;
            std::cin >> resp;
            if (resp == "e" || resp == "E") {
               std::cout << "expanding statistics at level "
                         << saved_sampler << std::endl;
               continue;
            }
            else if (resp == "r" || resp == "R") {
               std::cout << "reset statistics, repeat at level "
                         << saved_sampler << std::endl;
               sampler.reset_stats();
               continue;
            }
            else if (resp == "b" || resp == "B") {
               if (saved_sampler > 0) {
                  --saved_sampler;
               }
               std::stringstream sfile;
               sfile << "ex_" << saved_sampler << ".astate";
               sampler.restoreState(sfile.str().c_str());
               sampler.reset_stats();
               std::cout << "starting again at level "
                         << saved_sampler << std::endl;
               continue;
            }
         }

         std::stringstream sfile;
         sfile << "ex_" << saved_sampler << ".astate";
         sampler.saveState(sfile.str().c_str());
         old_efficiency = sampler.getEfficiency();
         if (sampler.adapt() > 0) {
            sampler.saveState("ex.astate");
            efficiency = sampler.getEfficiency();
            //sampler.display_tree();
            sampler.reset_stats();
         }
         if (old_efficiency > 0.1)
            sampler.setAdaptation_sampling_threshold(1000);
         std::cout << "advancing to level "
                   << ++saved_sampler << std::endl;
      }
      if (i == 500000000) {
         sum[0] = 0;
         sum[1] = 0;
         sum[2] = 0;
      }
   }
   sampler.saveState("ex.astate");
   sampler.display_tree();
   double error;
   double result = sampler.getResult(&error);
   std::cout << "result = " << result << " +/- " << error << std::endl;
   double mu = sum[1] / sum[0];
   double eff = sum[1] * sum[1] / (sum[0] * sum[2]);
   double rms = sqrt(sum[2] * (1 - eff)) / sum[0];
   std::cout << "result_IS is " << mu << " +/- " << rms 
             << ", efficiency = " << eff << std::endl;
   return 0;
}
