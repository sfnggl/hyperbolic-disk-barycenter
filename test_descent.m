addpath(genpath('./src'))

function test_descent
  max_iter = 10000;
  n = 5;
  %% Generate a positive-definite matrix
  A = rand(n,n);
  A = 0.5*(A+A');
  A += n*eye(n);
  %% Generate pseudo-random bias and stating point vectors
  b = rand(n)'; x0 = rand(n,1);
  %% Define function to evaluate
  func = @(x, A, b)(0.5*(x'*A*x) - b'*x);
  gradfunc = @(x, A, b)(A*x - b);
  %% easier call
  f = @(x)(func(x,A,b));
  g = @(x)(gradfunc(x,A,b));
  tic
  x_stationary = steepest_descent(f, g, n, max_iter, x0)
  toc
  A*x_stationary - b
end

