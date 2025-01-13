function d = partialder(x, a)
  % (R, R), (R, R) |-> (R, R)
  %
  % the very composite gradient of the distance function.
  %
  % define the single derivatives
  orange = 1 + 2 * (norm(x - a)^2 / (1-norm(x)^2) / (1-norm(a)^2));
  blue = 2*(x - a);
  if norm(x) == 0
    red = [0, 0]';
  else
    red = [x(1) / sqrt(x(1)^2 - x(2)^2), - x(2) / sqrt(x(1)^2 - x(2)^2)]';
  end
  pink = 2 ./ sqrt(orange.^2 - 1) .* ...
	 ((((1 - norm(x)^2) * (1 - norm(a)^2)) .* blue) + ...
	  ((norm(x - a)^2 * (1 - norm(a)^2) * 2 * norm(x)) .* red)) / ...
	 (((1 - norm(x)^2)*(1 - norm(a)^2))^2);
  d = pink;
end

function d = partial_derivative(z1, z2)
  orange = 1 + 2 * (norm(x - a)^2 / (1-norm(x)^2) / (1-norm(a)^2));
  
end

