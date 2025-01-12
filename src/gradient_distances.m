function G = gradient_distances(x, ps)
  G = 0;
  n = length(ps);
  for i = 1:n
    dis = distance(x,ps(:,i));
    de = partial_derivative(x,ps(:,i));
    G += dis*de*2/n;
  end
  G = G .* 2 ./ n;
end

function d = partial_derivative(z1, z2)
  g = @(z1, z2) (2 * norm(z1 - z2,1)^2);
  h = @(z1, z2) ((1 - norm(z1,1)^2)*(1 - norm(z2,1)^2));
  f = @(z1, z2) (1 + g(z1,z2) / h(z1,z2));
  dg = [4 * (z1 - z2), 4 * (z2 - z1)];
  dh = [z1 * (1 - norm(z2,1)^2), z2 * (1 - norm(z1,1)^2)] .* (-2);
  df = [...
    (1 - norm(z1,1)^2) * (1 - norm(z2,1)^2) * (z1 - z2) + z1 * norm(z1 - z2,1)^2 * (1 - norm(z2)^2), ... % d/d(z1)
    (1 - norm(z1,1)^2) * (1 - norm(z2,1)^2) * (z2 - z1) + z2 * norm(z2 - z1,1)^2 * (1 - norm(z1)^2), % d/d(z2)
  ];
  d = df * 4 / (sqrt(f(z1,z2)^2 - 1) * (1 - norm(z1,1)^2)^2 * (1 - norm(z2,1)^2)^2);
end
