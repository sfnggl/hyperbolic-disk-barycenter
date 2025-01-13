function b = poincarediskbaricenter()
% UNIPG ApNeA - A.A. 2024/25
%
% Method to search the baricenter of a Poincarè disk manifold
% It assumes n = 2
%
% it uses and compares
% -- steepest descent with Armijo's Method
% -- barzilai-borwein method with Newton's Method
% -- non-monotonic method with modified Newton's Method
% -- non-monotonic method with Wolfe conditions (TODO)
%
% This file is the entry point of the final project for the
% Numerical Approximation course at Unipg.
% Original authors: Stefano Gigli, Vito Festa, Giuliana Mammone, Carmine diMatteo

  clc;
  clear all;
  clf;
  pkg load sqlite;

  % read database
  db = sqlite('points.db');
  data_points = sqlread(db, "unit_circle_points");
  close(db);

  % initialize data structures
  points = [cell2mat(data_points.x),cell2mat(data_points.y)]';
  n = size(points,1);
  f = @(x) weighted_distances(x, points);
  g = @(x) gradient_distances(x, points);
  A = posdef(n); b = rand(n,1);
  % x0 = rand_start()
  x0 = [0;0]

  % initialize parameters
  if not(exist("max_iter"))
    max_iter = 100;
  end
  if not(exist("tol"))
    tol = 10d-10;
  end

  % plot of the points on the unit square
  arcs = linspace(0,2*pi,200)';
  circle_point = [cos(arcs), sin(arcs)]';

  % perform the searches
  xs = zeros(2*4, max_iter)';
  ds = zeros(2*4, max_iter)';
  ss = zeros(4, 1);
  bs = zeros(4, 1);
  steps = []; derivs = [];

  [steps, derivs, ss(1)] = sd(f,g,x0,max_iter,tol);
  xs([1:ss(1)],[1,2]) = steps';
  ds([1:ss(1)],[1,2]) = derivs';
  [steps, derivs, ss(2)] = bb(f,g,x0,max_iter,tol);
  xs([1:ss(2)],[3,4]) = steps';
  ds([1:ss(2)],[3,4]) = derivs';
  [steps, derivs, ss(3)] = nmt(f,g,x0,max_iter,tol);
  xs([1:ss(3)],[5,6]) = steps';
  ds([1:ss(3)],[5,6]) = derivs';
  [steps, derivs, ss(4)] = wolfe(f,g,x0,max_iter,tol);
  xs([1:ss(4)],[7,8]) = steps';
  ds([1:ss(4)],[7,8]) = derivs';

  % plot
  figure(1);
  hold on;
  ## subplot(1,4,1)
  plot([1:ss(1)], vecnorm(ds(:,[1,2]),2,2)(1:ss(1)), '-', 'LineWidth', 2, 'DisplayName', 'Steepest Descent', 'Color', 'green');
  legend();
  figure(2);
  ## subplot(1,4,2)
  plot([1:ss(2)], vecnorm(ds(:,[3,4]),2,2)(1:ss(2)), '-', 'LineWidth', 2, 'DisplayName', 'Barzilai-Borwein', 'Color', 'red');
  legend();
  figure(3);
  ## subplot(1,4,3)
  plot([1:ss(3)], vecnorm(ds(:,[5,6]),2,2)(1:ss(3)), '-', 'LineWidth', 2, 'DisplayName', 'Non Monotonic Search', 'Color', 'blue');
  legend();
  figure(4);
  ## subplot(1,4,4)
  plot([1:ss(4)], vecnorm(ds(:,[7,8]),2,2)(1:ss(4)), '-', 'LineWidth', 2, 'DisplayName', "Wolfe\'s Condition", 'Color', 'magenta');
  legend();
  hold off;

end

function [x,d,s] = dummy(n)
  x = randn(n,2);
  d = randn(n,2);
  s = n;
end

function x0 = rand_start()
  x = inf; y = inf;
  while (x^2 + y^2 >= 1)
    x = rand(); y = rand();
  end
  x0 = [x; y];
end

function D = weighted_distances(x, ps)
  % weighted sum of a fixed amount of points from the given x
  % C, (C, ... , C) |-> C
  %
  % it assumes equal mass for each point on the hyperbolic plane
  % argmin(n^-1 * sum(d(x,p_i)^2)) is solved for the baricenter x_b.
  %
  D = 0;
  n = length(ps);
  for i = 1:n
    D = D + distance(x,ps(:,i))^2;
  end
  D = D / n;
end

function G = gradient_distances(x, ps)
  % derivative of the weighted sums
  % C, (C, ... , C) |-> C
  %
  % it sums the partial derivative of the point
  % with respect of each point
  % it uses only the first component as
  % the direction for the next iteration
  % 
  % it's an highly unstable function for points
  % "a little too far" from the cluster of points
  G = 0;
  n = length(ps);
  for i = 1:n
    dis = distance(x,ps(:,i));
    de = partial_derivative(x,ps(:,i))(:,1);
    G += 2 / n^2 * dis * de;
  end
end

function d = partial_derivative(z1, z2)
  % the 2-dimensional gradient
  % with respect to the two points
  % C, C |-> C
  %
  g = @(z1, z2) (2 * norm(z1 - z2,1)^2);
  h = @(z1, z2) ((1 - norm(z1,1) ^2)  * (1 - norm(z2,1) ^2));
  f = @(z1, z2) (1 + g(z1,z2)  / h(z1,z2));
  dg = [4 * (z1 - z2), 4 * (z2 - z1)]; % (d/d(z1) g, d/d(z2) g)
  dh = [z1  * (1 - norm(z2,1) ^2) * (-2), z2  * (1 - norm(z1,1) ^2) * (-2)]; % (d/d(z1) h, d/d(z2) h)
  df = (h(z1,z2)  * dg - g(z1,z2)  * dh)  / h(z1,z2) ^2; % (d/d(z1) f, d/d(z2) f)
  d = df / sqrt(f(z1,z2)  ^ 2 - 1);
end

function d = distance(x, y)
  % the distance function on the hyperbolic plane
  % C, C |-> C
  %
  % derived from the equality which holds that
  % for each pair of points x,y in D, hyperbolic plane of n-dimension
  % norm(x-y)^2 / (1 - norm(x)^2) / (1 - norm(x)^2 =  5 * (cosh(d(x,y) - 1))
  %
  modsqrx = norm(x,1) ^ 2;
  modsqry = norm(y,1) ^ 2;
  esqrdiff = norm(x-y,1) ^ 2;
  d = acosh(1 + 2 * esqrdiff / (1 - modsqrx) / (1 - modsqry));
end

function [xs, ds, steps] = sd(f,g,x0, max_iter, tol)
  % Steepest Descent with Armijo's Rule for inexact line search
  % It uses a small starting α to offset the rapidly varying derivative
  % Iteratively increases the step taken to update the x with the gradient
  %
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [armijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(d(:,l)) > tol && l < max_iter
    % compute gradient
    x_next = x(:,l) + alpha(:,l)*d(:,l);
    if (norm(x_next,1) >= 1)
      break
    end
    x = [x, x_next];
    d = [d,-g(x(:,l+1))];
    % find new alpha via Armijo's method
    % by which chosing α s.t. max[f(x + α*d) < f(x) - σ*d'*d*α]
    % achieves convergence to the result
    if norm(d(:,l+1) - d(:,l)) <= tol
      break
    end
    alpha = [alpha,armijo(f,g,x(:,l+1))];
    l = l+1;
  end
  steps = l;
  xs = x;
  ds = d;
end

function alpha = armijo(f,g,x,opts)
  % backtracking line search
  %
  % given an adequate step size α, it can be proved
  % that a search over the condition
  % f(x + α*p) < f(x) + α*τ
  % given f ∈ C¹,τ = c*m, c ∈ (0,1), m > 0
  % yields a result x' in a finite amount of steps
  % such that ▽f(x') = 0
  %
  %
  if not(exist("opts"))
    abar = 0.1;
    beta = 0.4;
    sigma = 0.9;
  else
    abar = opts.a;
    beta = opts.b;
    sigma = opts.s;
  end
  alpha = abar;
  m = 0;
  d = -g(x);
  while f(x + alpha*d) >= f(x) - sigma*d'*d*alpha
    m++;
    if m > 500
      break
    end
    alpha = alpha * beta;
  end
end

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

function [xs, ds, steps] = wolfe(f,g,x0, max_iter, tol)
  % Steepest Descent with modified Armijo's Rule for inexact line search
  %
  % Iteratively reduces the step taken to update the x with the gradient
  %
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [modiarmijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(d(:,l)) > tol && l < max_iter
    x_next = x(:,l) + alpha(:,l)*d(:,l);
    if (norm(x_next,1) >= 1)
      break
    end
    x = [x, x_next];
    d = [d,-g(x(:,l+1))];
    if norm(d(:,l+1) - d(:,l)) <= tol
      break
    end
    % find new alpha via modified Armijo's condition
    % such that α
    % achieves convergence to the result
    alpha = [alpha,modiarmijo(f,g,x(:,l+1))];
    l = l+1;
  end
  steps = l;
  xs = x;
  ds = d;
end

function alpha = modiarmijo(f,g,x,opts)
  % modified armijo's condition
  %
  % using two bounds 0 < c1 < c2 < 1, an adequate
  % step size α which holds:
  % i) f(x - α*g) ≤ f(k) -c1*α*g'*f'(x)
  % ii) -g'*f'(x - α*g) ≤ c2*g'*f'(x)
  % reduces the objective function f
  %
  if not(exist("opts"))
    % Nocedal, Jorge; Wright, Stephen (1999). Numerical Optimization. p. 38.
    c1 = 10d-4;
    c2 = 0.9;
  else
    c1 = opts.c1;
    c2 = opts.c2;
  end
  assert(0 < c1 && c1 < c2 && c2 < 1);
  alpha = .1;
  p = -g(x);
  while f(x + alpha*p) > f(x) + c1*alpha*p'*g(x) && -p'*g(x + alpha*p) <= -c2*p'*g(x)
    alpha++;
    if alpha > 100
      break
    end
  end
end
