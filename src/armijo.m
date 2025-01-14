function alpha = armijo(f,g,x,opts)
  % backtracking line search
  %
  % given an adequate step size α, it can be proved
  % that a search over the condition
  % f(x + α*p) < f(x) + α*τ
  % given f ∈ C¹,τ = c*m, c ∈ (0,1), m > 0
  % yields a result x' in a finite amount of steps
  % such that ▽f(x') = 0
  %
  %
  if not(exist("opts"))
    abar = 0.1;
    beta = 0.4;
    sigma = 0.9;
  else
    abar = opts.a;
    beta = opts.b;
    sigma = opts.s;
  end
  alpha = abar;
  m = 0;
  d = -g(x);
  while f(x + alpha*d) >= f(x) - sigma*d'*d*alpha
    m++;
    if m > 500
      break
    end
    alpha = alpha * beta;
  end
end
