function b = poincare(n)
  %% random points on the unit square
  %% the evaluation function and its gradient
  %% perform the search
  %% plot results
end

function d = distance(p0, p1)
  %% the distance function
  d = acosh(1 + 2*pow(norm(P1 - p2),2) / (pow(1 - norm(P1),2) + pow(1 - norm(P2),2)))
end

function b = baricenter(x, w, p)
  %% the baricenter function
  b = 0
  for i = 1:size(w,1)
    b += w(i)*distance(x,p(i))
  end
end
