#ifndef hawkes_neg_hpp
#define hawkes_neg_hpp
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/* Adapted from Stelfi */
/* Marked Univariate Hawkes Process with Negative Alpha */
template<class Type>
Type hawkes_neg(objective_function<Type>* obj) {
  using namespace Eigen;
  
  // Data
  DATA_VECTOR(times);        
  DATA_VECTOR(marks);        
  DATA_SCALAR(penalty_coef); 
  int T = times.size();
  Type mean_marks = marks.sum()/marks.size();
  
  // Parameters
  PARAMETER(log_mu);         
  PARAMETER(a_par);     // for negative alpha model
  PARAMETER(log_beta);       
  
  // Transform parameters
  Type mu = exp(log_mu);
  Type beta = exp(log_beta);

  // Initialize A vector A(1)=0
  vector<Type> A = vector<Type>::Zero(T);
  // Calculate A vector with marks
  for(int j = 1; j < T; ++j) {
    A(j) = exp(-beta * (times[j] - times[j-1])) * (marks[j-1] + A(j-1));
    // A(i) = e^{-\beta(t_i-t_{i-1})}(M(i-1) + A(i-1))
  }
  
  // Calculate B vector
  vector<Type> B = vector<Type>::Zero(T);
  for(int i = 0; i < T; ++i) {
    B[i] = A[i] + marks[i];
  }
  
  // Calculate alpha with mark scaling
  // Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * (beta/mean_marks); 
  // enforcing 0<=alpha<=beta
  // -min(lambda/B) <= alpha <= beta/mean(marks)
  Type Max = beta/mean_marks;
  Type Min = mu/max(B);
  Type alpha = Type(0.5) * (Max + Min) * (tanh(a_par) + Type(1.)) - Min;

  // Calculate negative log-likelihood
  Type nll = Type(0.);
  Type last = times[T-1]; // t_n
  
  // Efficient marks sum calculation
  vector<Type> term_3vec = log(mu + alpha * A);
  nll = (mu * last)  // Baseline compensator
      + ((alpha / beta) * (marks.sum() - marks[T-1] - A(T-1)))  // Triggering compensator, Ozaki's efficient compensator algorithms
      - sum(term_3vec); // λ(t_i) = μ + α∑_{t_j<t_i}m(t_j)e^{-β(t_i-t_j)}
  
  // Regularization term on branching ratio

    nll += penalty_coef * (alpha / beta);

  
  ADREPORT(mu);
  ADREPORT(alpha);
  ADREPORT(beta);
  
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
