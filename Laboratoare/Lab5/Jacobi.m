function [x p] = Jacobi(A, b, x0, tol, max_iter)
  % diag(matrice) = vector cu elem de pe diagonala matricei
  % diag(vector) = matrice diagonala care are pe diagonala elementele din vector
  D = diag(diag(A));
  L_U = A - D;  % L + U

  % -D^-1 * (L + U)
  G = 1./D * L_U; % ASA NU din cauza 1/0
  G = diag(1./diag(D)) * L_U;

  % eig(G) = vector cu valorile proprii ale lui G
  RoG = max(abs(eig(G)));  % raza spectrala a lui G
  
  if RoG > 1
    x = NaN;
    step = -1;
    return;
  endif
  x = x0;

  for p = 1 : max_iter
    % x^(p) = x0
    % x^(p + 1) = D^-1 * (b - (L + U) * x^(p))
    x = diag(1./diag(D)) * (b - L_U * x0);  % x^(p + 1) = x

    if norm(x - x0) < tol
      break;
    endif

    x0 = x;
  endfor
endfunction