function A = posdef(n)
  a = 10^0
  b = 10^3
  D = diag(linspace(a,b,n));
  [Q,~] = qr(randn(n));
  A = Q * D * Q';
end

