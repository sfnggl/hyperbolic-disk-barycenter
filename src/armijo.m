function alpha = armijo(f,g,x,opts)
  if not(exist("opts"))
    abar = 1;
    beta = 0.5;
    sigma = 0.5;
  else
    abar = opts.a;
    beta = opts.b;
    sigma = opts.s;
  end
  alpha = abar;
  m = 1;
  d = -g(x);
  while f(x + alpha*d) >= f(x) - sigma*alpha*d'*d
    m++;
    alpha = alpha * beta;
  end
end

%% Si puo creare un file unico dove si
%% passano una struttura dati di opzioni
%% e in base a quelli ritornare una formula diversa
