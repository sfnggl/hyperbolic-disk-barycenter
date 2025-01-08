function alpha = bb(f, g, x0, max_iter, tol)
  if not(exist("max_iter"))
    max_iter = 1000
  end
  if not(exist("tol"))
    tol = 10d-10
  end
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [1];
  while norm(d(:,l)) > tol && l < max_iter
    x = [x, x + alpha(l)*d(:,l)];
    d = [d, -g(x(:,l+1))];
    s = x(:,l+1) - x(:,l);
    y = d(:,l+1) - d(:,l);
    alpha = [alpha, (s' * s) / (s' * y)];
    if norm(d(:,l) - d(:,l-1)) < tol
      break
    end
    l++;
  end
end
