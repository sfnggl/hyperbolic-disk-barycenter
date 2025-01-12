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
