function [x step] = SOR(A, b, x0, w, tol, max_iter)
  D = diag(diag(A));
  L = tril(A) - D;
  U = triu(A) - D;
  
  D_L = (D + w * L);
  U_D = (w * U + (w - 1) * D);

  G = D_L^-1 * U_D;

  if max(abs(eig(G))) > 1 + eps
    x = NaN;
    step = -1;
    return;
  endif

  n = length(b);
  x = x0;

  for step = 1 : max_iter
    % x = (D + L)^-1 * (b - (L + U) * x0)
    x = D_L^-1 * (b * w - U_D * x0);

    if norm(x - x0) < tol
      break;
    endif

    x0 = x;
  endfor
endfunction