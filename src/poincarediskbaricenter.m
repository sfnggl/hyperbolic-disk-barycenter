function b = poincarediskbaricenter(xs)
% UNIPG ApNeA - A.A. 2024/25
%
% Method to search the baricenter of a Poincarè disk manifold
% It assumes n = 2
%
% it uses and compares
% -- steepest descent with Armijo's Method
% -- barzilai-borwein method with Newton's Method
% -- non-monotonic method with modified Newton's Method
% -- non-monotonic method with Wolfe conditions (TODO)
%
% This file is the entry point of the final project for the
% Numerical Approximation course at Unipg.
% Original authors: Stefano Gigli, Vito Festa, Giuliana Mammone, Carmine diMatteo

  points = xs;

  % initialize data structures
  n = size(points,1);
  f = @(x) weighted_distances(x, points);
  g = @(x) gradient_distances(x, points);
  A = posdef(n);

  x0 = [0;0]

  % initialize parameters
  if not(exist("max_iter"))
    max_iter = 500;
  end
  if not(exist("tol"))
    tol = 10d-10;
  end

  % perform the searches
  xs = zeros(2*4, max_iter)';
  ds = zeros(2*4, max_iter)';
  ss = zeros(4, 1);
  bs = zeros(4, 2);
  steps = []; derivs = [];

  [steps, derivs, ss(1)] = gd(f,g,x0,max_iter,tol);
  xs([1:ss(1)],[1,2]) = steps';
  ds([1:ss(1)],[1,2]) = derivs';
  bs(1,:) = xs(ss(1),[1:2]);
  [steps, derivs, ss(2)] = bb(f,g,x0,max_iter,tol);
  xs([1:ss(2)],[3,4]) = steps';
  ds([1:ss(2)],[3,4]) = derivs';
  bs(2,:) = xs(ss(2),[3:4]);
  [steps, derivs, ss(3)] = nmt(f,g,x0,max_iter,tol);
  xs([1:ss(3)],[5,6]) = steps';
  ds([1:ss(3)],[5,6]) = derivs';
  bs(3,:) = xs(ss(3),[5:6]);
  [steps, derivs, ss(4)] = wolfe(f,g,x0,max_iter,tol);
  xs([1:ss(4)],[7,8]) = steps';
  ds([1:ss(4)],[7,8]) = derivs';
  bs(4,:) = xs(ss(4),[7:8]);

  % plot
  figure(1);
  hold on;
  ## subplot(1,4,1)
  semilogy([1:ss(1)], vecnorm(ds(:,[1,2]),2,2)(1:ss(1)), '-', 'LineWidth', 2, 'DisplayName', 'Gradient Descent', 'Color', 'green');
  title('Gradient Descent', 'FontSize', 18);
  xlabel('Steps', 'FontSize', 15);
  ylabel('||G||', 'FontSize', 15);
  set(gca, 'FontSize', 12);
  legend();
  figure(2);
  ## subplot(1,4,2)
  semilogy([1:ss(2)], vecnorm(ds(:,[3,4]),2,2)(1:ss(2)), '-', 'LineWidth', 2, 'DisplayName', 'Barzilai-Borwein', 'Color', 'red');
  title('Barzilai-Borwein', 'FontSize', 18);
  xlabel('Steps', 'FontSize', 15);
  ylabel('||G||', 'FontSize', 15);
  set(gca, 'FontSize', 12);
  legend();
  figure(3);
  ## subplot(1,4,3)
  semilogy([1:ss(3)], vecnorm(ds(:,[5,6]),2,2)(1:ss(3)), '-', 'LineWidth', 2, 'DisplayName', 'Non Monotonic Search', 'Color', 'blue');
  title('Non-Monotonic Search', 'FontSize', 18);
  xlabel('Steps', 'FontSize', 15);
  ylabel('||G||', 'FontSize', 15);
  set(gca, 'FontSize', 12);
  legend();
  figure(4);
  ## subplot(1,4,4)
  semilogy([1:ss(4)], vecnorm(ds(:,[7,8]),2,2)(1:ss(4)), '-', 'LineWidth', 2, 'DisplayName', "Wolfe\'s Conditions", 'Color', 'magenta');
  title("Wolfe\'s Conditions", 'FontSize', 18);
  xlabel('Steps', 'FontSize', 15);
  ylabel('||G||', 'FontSize', 15);
  set(gca, 'FontSize', 12);
  legend();
  hold off;

  b = bs;

end
