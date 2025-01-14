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
