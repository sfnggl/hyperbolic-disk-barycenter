function b = poincarediskbaricenter(opts)
% UNIPG ApNeA - A.A. 2024/25
%
% Method to search the baricenter of a Poincarè disk manifold
%
% It assumes n = 2
%
% depending on which optimization method is preferred:
% either steepest descent by means of Armijo's Method or
% barzilai-borwein method
% the second structure may define some hyperparameters used by Armijo's method
% though standard values (1, 0.5, 0.5) are ofter preferred.
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
  
  % perform the searches
  % sd
  % bb with alternating steps
  % bb with non-monotone line search
  [xs, ds, steps] = search(f,g,x0,opts.search_type);
  b = xs(:,steps)
  % plot results
  % plot points on the unit square (first subplot)
  % plot(unit circle, '-');
  % plot(points, 0, '*');
  % plot(xs, 0, '+');
  % plot convergence of the various methods (second subplot)
  semilogy([1:steps], ds, '-');

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
  D = 0;
  for i = 1:length(ps)
    D = D + distance(x,ps(i))^2;
  end
  D = D / n;
end

function G = gradient_distances(x, ps)

end

function d = distance(x, y)
  %% the distance function
  modsqrx = norm(x)^2;
  modsqry = norm(y)^2;
  esqrdiff = norm(x-y)^2;
  d = acosh(1 + 2 * esqrdiff / (1 - modsqrx) / (1 - modsqry));
end
