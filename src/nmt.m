function [xs, ds, steps] = nmt(f,g,x0, max_iter, tol)
  % Non-Monotonic line search
  %
  % usage of modified Netwon condition
  % higher amount of backtracking permitted
  % provided it is bounded by a at most j-th previous value
  % scaled by some scalar ∈ (0,1) and
  % traslated by a constant C,
  % linearly dependent on a parameter η and the value of the next function
  %
  M = 10d3;
  gamma = .5;
  sigma = .9;
  alpha = .1;
  theta = 1;
  l = 1;
  k = 1;
  x = [x0];
  d = [-g(x0)];
  t = [f(x0)];
  m = [1];
  while norm(d(:,l)) > tol && l < max_iter
    x_next = x(:,l) + alpha*d(:,l);
    if (norm(x_next,1) >= 1)
      break
    end
    x = [x, x_next];
    d = [d,-g(x(:,l+1))];
    t = [t, f(x(:,l+1))];
    if norm(d(:,l+1) - d(:,l)) < tol
      l = l+1;
      break
    end
    % find new α by satisfying modified Newton's method
    % if such a step can no longer be decided scale α by σ
    newton_cond = 0;
    for i = 1:m(k)
      newton_cond += (t(:,l+1) <= t(:,i) - gamma * alpha * d(:,l)' * d(:,l));
    end
    if (newton_cond)
      k = k + 1;
      m(k) = min(m(k-1)+1, M);
    else
      alpha = alpha * sigma;
    end
    l = l + 1;
  end
  steps = l;
  xs = x;
  ds = d;
end
