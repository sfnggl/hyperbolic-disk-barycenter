function G = gradient_distances(x, ps)
  % derivative of the weighted sums
  % C, (C, ... , C) |-> C
  %
  % it sums the partial derivative of the point
  % with respect of each point
  % it uses only the first component as
  % the direction for the next iteration
  %
  G = 0;
  n = length(ps);
  for i = 1:n
    dis = distance(x,ps(:,i));
    de = partial_derivative2(x,ps(:,i))(:,1);
    G += 2 / n .* dis .* de;
  end
end
