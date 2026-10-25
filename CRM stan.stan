//
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// removed truncation for ID70 prior, and set SD of ID70 to 1.5

data {
  int<lower=0> K; // number of doses/dose levels
  vector[K] dose; // doses/dose levels
  int<lower=0> n[K]; // number of participants at each dose
  int<lower=0> y[K]; // number of participants infected
}

// The parameters accepted by the model. Our model
// accepts two parameters 'a' and 'b'.
parameters {
  // log10(ID70) ~ N(4, 1.5^2)
  real <lower=1, upper=9> log_id70;
  real <lower=log(0.3), upper=log(0.8)> log_a;
}

transformed parameters {
  real a;
  real b;
  real id70;
  // convert to original scale
  a = exp(log_a);
  
  id70 = pow(10, log_id70);
  
  b = id70 / (pow(0.3, -1/a) - 1);
    
  vector[K] p;
  p = 1 - pow(1 + dose / b, -a);
}



// likelihood and prior 
model {
  // priors
  log_id70 ~ normal(4, 0.8);
  
  // likelihood
  y ~ binomial(n, p);
}


