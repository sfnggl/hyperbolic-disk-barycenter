function b = poincare(opts)


  pkg load sqlite;
  if ~exist("opts")
    search_type = struct("search_type", "sd");
  end
  %% random points on the unit square
  db = sqlite('points.db');
  data_points = sqlread(db, "unit_circle_points");
  points = [cell2mat(data_points.x), cell2mat(data_points.y)]';
  n = size(points,1);
  %% initialize data structures
  ws = ones(n,1) ./ n;
  x0 = rand_start()
  %% the distance function and its gradient
  distfunc =
  graddistfunc =
  f = @(p) func
  g = @(p) gradfunc
  %% perform the search
  [xs, ds, steps] = search(f,g,x0,search_type);
  b = xs(:,steps)
  %% plot results
  semilogy([1:steps], ds, '-');

end

function x0 = rand_start()
  x = inf; y = inf;
  while (x^2 + y^2 >= 1)
    x = rand(); y = rand();
  end
  x0 = [x; y];
end

function d = distance(x, y)
  %% the distance function
  modsqrx = abs(x)^2;
  modsqry = abs(y)^2;
  esqrdiff = abs(x-y)^2;
  d = acosh(1 + 2 * esqrdiff / (1 - modsqrx) / (1 - modsqry));
end

function b = baricenter(x, ws, ps)
  %% the baricenter function
  b = 0;
  for i = 1:size(ws,1)
    b += ws(i)*(distance(x,ps(:,i))^2);
  end
end

function b = baricenterreal(x)
  cidizeta = @(x) ((1 - x)/(1 + z));
  b = nthroot(cidizeta(x),size(x,1));
end
