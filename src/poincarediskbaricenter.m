
function b = poincarediskbaricenter(opts)
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

  pkg load sqlite;

  if ~exist("opts","var")
    error("Invalid call to poincarediskbaricenter. Correct usage is:\n\n\
-- B = poincarediskbaricenter(struct(\"search_type\", \"...\"))\n\n\n");
    return;
  end

  % random points on the unit square
  db = sqlite('src/points.db');
  data_points = sqlread(db, "unit_circle_points");
  close(db);

  %% initialize data structures
  points = [cell2mat(data_points.x),cell2mat(data_points.y)]';
  n = size(points,1);
  f = @(x) weighted_distances(x, points);
  g = @(x) gradient_distances(x, points);
  x0 = rand_start()
  if not(exist("max_iter"))
    max_iter = 1000;
  end
  if not(exist("tol"))
    tol = 10d-10;
  end
  
  % perform the searches
  xs = []; ds = []; steps = [];
  % sd w/ Armijo
  [xs(1), ds(1), steps(1)] = sd(f,g,x0, max_iter, tol);
  % bb with alternating steps
  [xs(2), ds(2), steps(2)] = bb(f,g,x0, max_iter, tol);
  % non-monotone line search with newton's
  [xs(3), ds(3), steps(3)] = nmt(f,g,x0, max_iter, tol);
  % non-monotone line search with wolfe's
  ## [xs(4), ds(4), steps(4)] = wolfe(f,g,x0, max_iter, tol);
  
  % collect the results
  bs = xs(:,steps)

  % plot results
  arcs = linspace(0,2*pi,200)';
  circle_point = [cos(arcs), sin(arcs)];
  % plot of the points on the unit square
  figure(1);
  plot(circle_point(:,1), circle_point(:,2), '-');
  plot(linspace(-1,1,100),zeros(100));
  plot(points(:,1), points(:,2), '+', 'markersize', 5);
  plot(xs(:,1)(1), xs(:,1)(2), '--', 'markersize', 2, 'green');
  plot(xs(:,2)(1), xs(:,2)(2), '--', 'markersize', 2, 'red');
  plot(xs(:,3)(1), xs(:,3)(2), '--', 'markersize', 2, 'blue');
  ## plot(xs(:,2)(1), xs(:,2)(2), '--', 'markersize', 2, 'orange');
  % plot convergence of the various methods (second subplot)
  figure(2);
  semilogy([1:steps], ds(:,1), '-', 'green');
  semilogy([1:steps], ds(:,2), '-', 'red');
  semilogy([1:steps], ds(:,3), '-', 'blue');
  ## semilogy([1:steps], ds(:,4), '-', 'orange');
end

function x0 = rand_start(n)
  x = inf; y = inf;
  while (x^2 + y^2 >= 1)
    x = rand(); y = rand();
  end
  switch (n)
    case 1
      x0 = x;
    case 2
      x0 = [x, y];
  end
end

function D = weighted_distances(x, ps)
  % weighted sum of a fixed amount of points from the given x
  % C, (C, ... , C) |-> C
  % 
  % it assumes equal mass for each point on the hyperbolic plane
  % argmin(n^-1 * sum(d(x,p_i)^2)) is solved at the baricenter B.
  % 
  D = 0;
  n = length(ps);
  for i = 1:n
    D = D + distance(x,ps(i))^2;
  end
  D = D / n;
end

function d = distance(x, y)
  % the distance function on the hyperbolic plane
  % C, C |-> C
  %
  % derived from the equality which holds that
  % for each pair of points x,y in D, hyperbolic plane of n-dimension
  % norm(x-y)^2 / (1 - norm(x)^2) / (1 - norm(x)^2 = .5 * (cosh(d(x,y) - 1))
  %
  modsqrx = norm(x)^2;
  modsqry = norm(y)^2;
  esqrdiff = norm(x-y)^2;
  d = acosh(1 + 2 * esqrdiff / (1 - modsqrx) / (1 - modsqry));
end

function G = gradient_distances(x, ps)
  G = 0;
  n = length(ps);
  de = @(x, p) partgdis(x(1), x(2), p(1), p(2));
  for i = 1:n
    G = G + 2*distance(x,ps(i))*de(x,ps(i));
  end
  G = G ./ n;
end

function g = partgdis(x,y,p,q)
  % (R, R), (R, R) |-> (R, R)
  %
  % the very composite gradient of the distance function.
  %
  % define the single derivatives
  arccosarg = 1 + 2 * norm([x, y] - [p, q])^2 / ...
		  (1 - norm([x, y])^2) / ...
		  (1 - norm([p, q])^2);  
  dnormdiffx = (x - p) / sqrt((x - p)^2+(y - q)^2);
  dnormdiffy = (y - q) / sqrt((x - p)^2+(y - q)^2);
  dnormx = x / sqrt(x^2 - y^2);
  dnormy = y / sqrt(x^2 - y^2);
  % easier aliases
  orange = 1 / sqrt(arccosarg^2 - 1);
  blue = [dnormdiffx, dnormdiffy];
  pink = [dnormx, dnormy];
  n1 = norm([x,y]);
  n2 = norm([p,q]);
  middle = (1-n1^2) * (1-n2^2) - ...
	   norm([x,y] + 2*[p,q])^2 * (1-n2^2) * n1;
  bottom = ((1-n1^2)*(1-n2^2))^2;
  % compose (d/dx, d/dy)
  g = (2 * orange) * bottom * middle .* blue .* pink;
end

function [xs, ds, steps] = sd(f,g,x0, max_iter, tol)
  % Steepest Descent with Armijo's Rule for inexact line search
  %
  % Iteratively reduces the step taken to update the x with the gradient
  %
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [armijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(x(:,l)) < 1 && norm(d(:,l)) > tol && l < max_iter
    % compute gradient
    x = [x, x(:,l) + alpha(:,l)*d(:,l)];
    d = [d,-g(x(:,l+1))];
    % find new alpha via Armijo's method
    % by which chosing α s.t. max[f(x + α*d) < f(x) - σ*d'*d*α]
    % achieves convergence to the result
    alpha = [alpha,armijo(f,g,x(:,l+1))];
    l = l+1;
    if norm(d(:,l) - d(:,l-1)) < tol
      break
    end
  end
  steps = l;
  x_stat = x;
  g_stat = d;
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
    abar = 1;
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
  % alternates between to stable variations
  % of the step size to take, uses second-order information
  % by means of Newton's Method, whereby iterating
  % x(k+1) = x(k) − (F(k))^−1*g(k)
  % founds a rapid solution
  %
  l = 1;
  x = [x0];
  d = [-g(x0)];
  alpha = [1];
  while norm(d(:,l)) < 1 && norm(d(:,l)) > tol && l < max_iter
    x = [x, x(:,l) - alpha(l)*d(:,l)];
    d = [d, -g(x(:,l+1))];
    s = x(:,l+1) - x(:,l);
    y = d(:,l+1) - d(:,l);
    %% compute next α
    if mod(l, 2)
      alpha = [alpha, (s' * s) / (s' * y)];
    else
      alpha = [alpha, (s' * y) / (y' * y)];
    end
    if norm(d(:,l) - d(:,l+1)) < tol
      break
    end
    l = l+1;
  end
  steps = l;
  x_stat = x;
  g_stat = d;
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
  sigma = .9
  alpha = 1;
  theta = 1;

  l = 1;
  x = [x0];
  d = [-g(x0)];
  f = [f(x0)];
  m = [0];

  while norm(d(:,l)) < 1 && norm(d(:,l)) > tol && l < max_iter
    x = [x, x(:,l) + alpha*d(:,l)];
    d = [d, -g(x(:,l+1))];
    f = [f, f(x(:,l+1))];
    
    % find new α by satisfying modified Newton's method
    % if such a step can no longer be decided scale α by σ
    newton_cond = 0;
    for i = 1:m(:,l)
      newton_cond += (f(:,l+1) <= f(:,i) - gamma * alpha * g(:,l)' * g(:,l));
    end
    if (newton_cond)
      k = k + 1;
      m(k) = min(m(k-1)+1, M);
    else
      alpha = alpha * sigma;
    end
  end
end

function [xs, ds, steps] = wolfe(f,g,x0, max_iter, tol)
  % Steepest Descent with Armijo's Rule for inexact line search
  %
  % Iteratively reduces the step taken to update the x with the gradient
  %
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [modiarmijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(x(:,l)) < 1 && norm(d(:,l)) > tol && l < max_iter
    x = [x, x(:,l) + alpha(:,l)*d(:,l)];
    d = [d,-g(x(:,l+1))];
    % find new alpha via modified Armijo's condition
    % such that α = min[i) && ii)]. this α
    % achieves convergence to the result
    alpha = [alpha,modiarmijo(f,g,x(:,l+1))];
    l = l+1;
    if norm(d(:,l) - d(:,l-1)) < tol
      break
    end
  end
  steps = l;
  x_stat = x;
  g_stat = d;
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
  alpha = 1;
  p = -g(x);
  while f(x + alpha*p) > f(x) + c1*alpha*p'*g(x) && -p'*g(x + alpha*p) <= -c2*p'*g(x)
    alpha++;
    if m > 100
      break
    end
  end
end
