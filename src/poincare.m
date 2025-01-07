function b = poincare(n)
  f = @(x, A, b) 0.5*(x'*A*x) - b'*x;
  g = @(x, A, b) A*x - b;
  steepest_descent(f,g);
end
