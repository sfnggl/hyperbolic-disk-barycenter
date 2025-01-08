function alpha = armijo(f,g,x,opts)
  if not(exist("opts"))
    abar = 1;
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

%% Si puo creare un file unico dove si
%% passano una struttura dati di opzioni
%% e in base a quelli ritornare una formula diversa
