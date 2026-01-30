#ifndef hawkes_multi_neg_hpp
#define hawkes_multi_neg_hpp
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/* Multi-Series Marked Univariate Hawkes Process with Negative Alpha - Shared Parameters Only */
template<class Type>
Type hawkes_multi_neg(objective_function<Type>* obj) {
  using namespace Eigen;
  
  // Data
  DATA_VECTOR(times);          // Concatenated times across all series
  DATA_VECTOR(marks);          // Concatenated marks across all series  
  DATA_IVECTOR(series_id);     // Series identifier for each event (0-indexed)
  DATA_SCALAR(penalty_coef);   // Regularization coefficient
  DATA_INTEGER(n_series);      // Number of series
  
  // Parameters - always shared across all series
  PARAMETER(log_mu);           // Single mu for all series
  // PARAMETER(logit_abratio);    // Single alpha ratio for all series - for positive alpha
  PARAMETER(a_par);            // Single a_par for negative alpha for all series
  PARAMETER(log_beta);         // Single beta for all series
  
  Type total_nll = Type(0.);
  // Calculate overall mean marks across all series for shared parameters
  Type overall_mean_marks = marks.sum() / marks.size();
  
  // Transform shared parameters
  Type mu = exp(log_mu);
  Type beta = exp(log_beta);
  // Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * (beta / overall_mean_marks); // Commented for negative alpha
  // Calculate max(B) across all series for constraint
  Type overall_max_B = Type(0.);
  std::vector<vector<Type>> a_vectors(n_series);
  
  // First pass: Calculate max(B) and store A_s vectors
  for(int s = 0; s < n_series; s++) {
    // Collect indices for this series
    vector<int> indices_s;
    for(int i = 0; i < times.size(); i++) {
      if(series_id[i] == s) {
        indices_s.conservativeResize(indices_s.size() + 1);
        indices_s[indices_s.size() - 1] = i;
      }
    }
    
    int T_s = indices_s.size();
    if(T_s == 0) {
      a_vectors[s] = vector<Type>::Zero(0);
      continue; // Skip empty series
    }
    
    // Sort indices by time for this series
    for(int i = 0; i < T_s - 1; i++) {
      for(int j = i + 1; j < T_s; j++) {
        if(times[indices_s[i]] > times[indices_s[j]]) {
          int temp = indices_s[i];
          indices_s[i] = indices_s[j];
          indices_s[j] = temp;
        }
      }
    }
    
    // Initialize A vector for this series - A(1)=0
    vector<Type> A_s = vector<Type>::Zero(T_s);
    
    // Calculate A vector with marks for this series
    for(int j = 1; j < T_s; ++j) {
      int curr_idx = indices_s[j];
      int prev_idx = indices_s[j-1];
      A_s[j] = exp(-beta * (times[curr_idx] - times[prev_idx])) * (marks[prev_idx] + A_s[j-1]);
    }
    a_vectors[s] = A_s; // Store A_s
    
    // Calculate B vector for this series and update overall max
    for(int i = 0; i < T_s; ++i) {
      Type B_si = A_s[i] + marks[indices_s[i]];
      if(B_si > overall_max_B) {
        overall_max_B = B_si;
      }
    }
  }
  
  // Calculate shared alpha with constraint
  Type Max = beta/overall_mean_marks;
  Type Min = mu / (overall_max_B * Type(n_series));
  Type alpha = Type(0.5) * (Max + Min) * (tanh(a_par) + Type(1.)) - Min;
  
  // Second pass: likelihood calculations for each series
  for(int s = 0; s < n_series; s++) {
    // Collect indices for this series
    vector<int> indices_s;
    for(int i = 0; i < times.size(); i++) {
      if(series_id[i] == s) {
        indices_s.conservativeResize(indices_s.size() + 1);
        indices_s[indices_s.size() - 1] = i;
      }
    }
    
    int T_s = indices_s.size();
    if(T_s == 0) continue; // Skip empty series
    
    // Sort indices by time for this series
    for(int i = 0; i < T_s - 1; i++) {
      for(int j = i + 1; j < T_s; j++) {
        if(times[indices_s[i]] > times[indices_s[j]]) {
          int temp = indices_s[i];
          indices_s[i] = indices_s[j];
          indices_s[j] = temp;
        }
      }
    }
    
    // Retrieve the pre-calculated A_s vector
    vector<Type> A_s = a_vectors[s];
    
    // Calculate negative log-likelihood for this series
    int last_idx = indices_s[T_s-1];
    Type last_s = times[last_idx];
    
    // Series marks sum
    Type marks_sum_s = Type(0.);
    for(int j = 0; j < T_s; ++j) {
      marks_sum_s += marks[indices_s[j]];
    }
    
    vector<Type> term_3vec_s = log(mu + alpha * A_s);
    Type nll_s = (mu * last_s)  // Baseline compensator
        + ((alpha / beta) * (marks_sum_s - marks[last_idx] - A_s[T_s-1]))  // Triggering compensator, Ozaki's efficient compensator algorithms
        - sum(term_3vec_s); // λ(t_i) = μ + α∑_{t_j<t_i}m(t_j)e^{-β(t_i-t_j)}
    
    /* Original loop-based approach (commented out for efficiency)
    // Log intensity terms
    Type log_intensity_sum = Type(0.);
    for(int j = 0; j < T_s; ++j) {
      log_intensity_sum += log(mu + alpha * A_s[j]);
    }
    
    Type nll_s = (mu * last_s)  // Baseline compensator
      + ((alpha / beta) * (marks_sum_s - marks[last_idx] - A_s[T_s-1]))  // Triggering compensator
      - log_intensity_sum;
    */
    
    total_nll += nll_s;
  }
  
  // Add regularization on branching ratio (single penalty for shared alpha)

    total_nll += penalty_coef * (alpha / beta);

  
  // Report shared parameters
  ADREPORT(mu);
  ADREPORT(alpha);
  ADREPORT(beta);
  
  return total_nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif