function alpha = modiarmijo(f,g,x,opts)
  % modified armijo's condition
  %
  % using two bounds 0 < c1 < c2 < 1, an adequate
  % step size α which holds:
  % i) f(x - α*g) ≤ f(k) -c1*α*g'*f'(x)
  % ii) -g'*f'(x - α*g) ≤ c2*g'*f'(x)
  % reduces the objective function f
  %
  if not(exist("opts"))
    % Nocedal, Jorge; Wright, Stephen (1999). Numerical Optimization. p. 38.
    c1 = 10d-4;
    c2 = 0.9;
  else
    c1 = opts.c1;
    c2 = opts.c2;
  end
  assert(0 < c1 && c1 < c2 && c2 < 1);
  alpha = .1;
  p = -g(x);
  while f(x + alpha*p) > f(x) + c1*alpha*p'*g(x) && -p'*g(x + alpha*p) <= -c2*p'*g(x)
    alpha++;
    if alpha > 100
      break
    end
  end
end
