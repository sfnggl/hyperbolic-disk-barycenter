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

  if ~exist("opts")
    opts = struct("search_type", "sd", "field", "real");
  end

  % random points on the unit square
  db = sqlite('points.db');
  data_points = sqlread(db, "unit_circle_points");

  %% initialize data structures based on given domain
  switch (opts.real)
	 case {"real"}
	   % search the baricenter on the diameter of the disk, i.e. Im(z) == 0
	   points = [cell2mat(data_points.x)]';
	   n = size(points,1);
	   f = @(x) sum(cellfun(@(pi) distance(x,pi)^2, points)) / n;
	   g = @(x) sum(cellfun(@(pi) (distance(x,pi) / (x^2 - 1)), points)) ...
	       * 2 / n;
	   x0 = rand_start()
	 case {"complex"}
  	  %
	 otherwise
	   disp(...
		 "Please provide a correct domain.
		 Possible domains are: real, complex")
	   exit()
  end
  
  % perform the search
  [xs, ds, steps] = search(f,g,x0,opts.search_type);
  b = xs(:,steps)
  % plot results
  % plot(unit circle, '-');
  % plot(points, 0, '*');
  % plot(xs, 0, '+');
  semilogy([1:steps], ds, '-');

end

function x0 = rand_start(n)
  if ~exist("n", "var")
    n = 1;
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

function d = distance(x, y)
  %% the distance function
  modsqrx = norm(x)^2;
  modsqry = norm(y)^2;
  esqrdiff = norm(x-y)^2;
  d = acosh(1 + 2 * esqrdiff / (1 - modsqrx) / (1 - modsqry));
end
