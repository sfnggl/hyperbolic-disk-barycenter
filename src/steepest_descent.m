function x_stat = steepest_descent(f, g, n, max_iter, x)
  %% initial step
  d = -g(x);
  alpha = armijo(f,g,x);
  i = 0;
  %% Iterate up to a numerical bound or convergence
  while (d != 0 && i < max_iter)
    %% compute gradient
    x = x + alpha*d;
    d = -g(x);
    alpha = armijo(f,g,x);
    i += 1;
  end
  x_stat = x;
end
