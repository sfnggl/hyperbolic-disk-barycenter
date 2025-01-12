function g = derivatives(x,y,p,q)
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
  blue = [dnormdiffx; dnormdiffy];
  pink = [dnormx; dnormy];
  n1 = norm([x,y]);
  n2 = norm([p,q]);
  middle = (1-n1^2) * (1-n2^2) - ...
	   norm([x,y] + 2*[p,q])^2 * (1-n2^2) * n1;
  bottom = ((1-n1^2)*(1-n2^2))^2;
  % compose (d/dx, d/dy)
  g = (2 * orange) * bottom * middle .* blue .* pink;
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

