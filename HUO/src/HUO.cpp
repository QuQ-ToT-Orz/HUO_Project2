#include <fstream>
#include <TMB.hpp>
#include <vector>
#include <iostream>
#include <numeric>
#include <math.h>

#include "hawkes.h"
#include "hawkes_neg.h"
#include "hawkes_multi.h"
#include "hawkes_multi_neg.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model_type);
  if (model_type == "hawkes") {
    return hawkes(this);
  } else
  if (model_type == "hawkes_neg") {
    return hawkes_neg(this);
  } else
  if (model_type == "hawkes_multi") {
    return hawkes_multi(this);
  } else
  if (model_type == "hawkes_multi_neg") {
    return hawkes_multi_neg(this);
  } else {
	Rf_error("Unknown model.");
   }
   return 0;
}
