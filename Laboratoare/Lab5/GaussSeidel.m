function [x p] = GaussSeidel(A, b, x0, tol, max_iter)
  % triu(A) = D + U
  D_L = tril(A);  % D + L
  U = A - D_L;  % A - D -L

  % G = -(D + L)^-1 * U
  G = -D_L^-1 * U;
  
  % eig(G) = vector cu valorile proprii ale lui G
  RoG = max(abs(eig(G)));  %raza spectrala a lui G

  if RoG > 1
    x = NaN;
    step = -1;
    return;
  endif

  x = x0;

  for p = 1 : max_iter
    % x^(p) = x0
    % x = (D + L)^-1 * (b - U * x0)
    x = D_L^-1 * (b - U * x0);  % x^(p + 1) = x

    if norm(x - x0) < tol
      break;
    endif

    x0 = x;
  endfor
endfunction