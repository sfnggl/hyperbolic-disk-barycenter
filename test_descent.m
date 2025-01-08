clear all;

addpath(genpath('./src'))

function [A, x_stationary] = test_descent
  max_iter = 1000;
  n = 10;
  %% Generate a positive-definite matrix
  A = rand(n,n);
  A = 0.5*(A+A');
  A += n*eye(n);
  %% Generate pseudo-random bias and stating point vectors
  b = rand(n,1); x0 = rand(n,1);
  %% Define function to evaluate
  func = @(x, A, b)(0.5*(x'*A*x) - b'*x);
  gradfunc = @(x, A, b)(A*x - b);
  %% easier call
  f = @(x)(func(x,A,b));
  g = @(x)(gradfunc(x,A,b));
  tic
  [xs_stationary, grad_sequence] = steepest_descent(f, g, x0);
  toc
  total_steps = size(xs_stationary, 2);
  semilogy([1:total_steps], cellfun(@norm,num2cell(grad_sequence,1)), "+", "markersize", 20);
end

