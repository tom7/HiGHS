#ifndef __SRC_LIB_PRICING_STEEPESTEDGEPRICING_HPP__
#define __SRC_LIB_PRICING_STEEPESTEDGEPRICING_HPP__

#include "pricing.hpp"
#include "../runtime.hpp"
#include "../basis.hpp"

// 

class SteepestEdgePricing : public Pricing {
private:
   Runtime& runtime;
   Basis& basis; 

   std::vector<double> weights;

   int chooseconstrainttodrop(const Vector& lambda) {
      auto activeconstraintidx = basis.getactive();
      auto constraintindexinbasisfactor = basis.getindexinfactor(); 
      
      int minidx = -1;
      double maxval = 0.0;
      for (unsigned int i = 0; i < activeconstraintidx.size(); i++) {
         int indexinbasis = constraintindexinbasisfactor[activeconstraintidx[i]];
         if (indexinbasis == -1) {
            printf("error\n");
         }
         assert(indexinbasis != -1);

         double val = lambda.value[indexinbasis] * lambda.value[indexinbasis] / weights[indexinbasis];
         if (val > maxval && fabs(lambda.value[indexinbasis]) > runtime.settings.lambda_zero_threshold) {
            if (basis.getstatus(activeconstraintidx[i]) == BasisStatus::ActiveAtLower &&
               -lambda.value[indexinbasis] > 0) {
               minidx = activeconstraintidx[i];
               maxval = val;
            } else if (basis.getstatus(activeconstraintidx[i]) == BasisStatus::ActiveAtUpper &&
                        lambda.value[indexinbasis] > 0) {
               minidx = activeconstraintidx[i];
               maxval = val;
            } else {
               // TODO
            }
         }
         
      }

      return minidx;
   }

public:  
   SteepestEdgePricing(Runtime& rt, Basis& bas) : runtime(rt), basis(bas), weights(std::vector<double>(rt.instance.num_var, 1.0)) {};
   
   int price(const Vector& x, const Vector& gradient) {
      Vector lambda = basis.ftran(gradient);
		int minidx = chooseconstrainttodrop(lambda);
      return minidx;
   }

   void update_weights(const Vector& aq, const Vector& ep, unsigned int p, unsigned int q) {
      unsigned int rowindex_p = basis.getindexinfactor()[p];
      
      Vector v = basis.btran(aq);
      
      double weight_p = weights[rowindex_p];
      for (unsigned int i=0; i<runtime.instance.num_var; i++) {
         if (i == rowindex_p) {
            weights[i] = weight_p / (aq.value[rowindex_p] * aq.value[rowindex_p]);
         } else {
            weights[i] = weights[i] - 2 * (aq.value[i] / aq.value[rowindex_p]) * (v.value[i])   + (aq.value[i] * aq.value[i]) / (aq.value[rowindex_p] * aq.value[rowindex_p]) * weight_p;
         }
      }
   }

};

#endif
