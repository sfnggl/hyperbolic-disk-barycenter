function [xs, ds, steps] = bb(f,g,x0, max_iter, tol)
  % Barzilai-Borwein with Big-Small step
  %
  % uses second-order information
  % by means of Newton's Method, whereby iterating
  % x(k+1) = x(k) − (F(k))^−1*g(k)
  % founds a rapid solution
  %
  l = 1;
  x = [x0];
  d = [-g(x0)];
  alpha = [0.1];
  while norm(d(:,l)) > tol && l < max_iter
    x_next = x(:,l) + alpha(:,l)*d(:,l);
    if (norm(x_next,1) >= 1)
      break
    end
    x = [x, x_next];
    d = [d,-g(x(:,l+1))];
    s = x(:,l+1) - x(:,l);
    y = d(:,l+1) - d(:,l);
    if norm(d(:,l+1) - d(:,l)) <= tol
      break
    end
    %% compute next α
    alpha = [alpha, -1 * (s' * s) / (s' * y)];
    l = l+1;
  end
  steps = l;
  xs = x;
  ds = d;
end
