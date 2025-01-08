function b = poincare(n,opts)
  if ~exist("n") || n <= 0
    error("no amount of points supplied")
  end
  if ~exist("opts")
    search_type = opts.search_type
  else
    search_type = "sd"
  end
  %% random points on the unit square
  db = sqlite('points.db')
  data_points = sqlread(db, "unit_circle_points")
  points = table2array(data_points)
  %% the evaluation function and its gradient
  func = @(x,y) x^2 + y^2
  gradfunc = @(x,y) 2*x + 2*y
  f = @(p) func(p(0),p(1))
  g = @(p) gradfunc(p(0),p(1))
  %% initialize data structures
  bs = [];
  ws = ones(n,1)./n;
  x0 = rand_start()
  %% perform the search
  [xs, ds, steps] = search(f,g,x0,search_type);
  b = xs(steps)
  %% plot results
  plot([1:steps], baricenter(xs, ws, points), '-');
end

function x0 = rand_start()
  x = inf; y = inf
  while (x^2 + y^2 >= 1)
    x = rand(); y = rand()
  end
  x0 = [x, y]
end

function d = distance(p0, p1)
  %% the distance function
  d = acosh(1 + 2*pow(norm(P1 - p2),2) / (pow(1 - norm(P1),2) + pow(1 - norm(P2),2)))
end

function b = baricenter(xs, ws, ps)
  %% the baricenter function
  b = 0
  for i = 1:size(ws,1)
    b += ws(i)*distance(xs,ps(i))
  end
end

function b = baricenterreal(xs)
  cidizeta = @(xs) (1-xs)./(1+z)
  b = nthroot(cidizeta(xs),size(xs,1))
end
