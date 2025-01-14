function d = partial_derivative(z1, z2)
  % the 2-dimensional gradient
  % with respect to the two points
  % C, C |-> C
  %
  g = @(z1, z2) (2 * norm(z1 - z2,1)^2);
  h = @(z1, z2) ((1 - norm(z1,1) ^2)  * (1 - norm(z2,1) ^2));
  f = @(z1, z2) (1 + g(z1,z2) / h(z1,z2));
  dg = [4 * (z1 - z2), 4 * (z2 - z1)]; % (d/d(z1) g, d/d(z2) g)
  dh = [z1  * (1 - norm(z2,1) ^2) * (-2), z2  * (1 - norm(z1,1) ^2) * (-2)]; % (d/d(z1) h, d/d(z2) h)
  df = (h(z1,z2) * dg - g(z1,z2) * dh)  / h(z1,z2) ^2; % (d/d(z1) f, d/d(z2) f)
  d = df / sqrt(f(z1,z2) ^ 2 - 1);
end
