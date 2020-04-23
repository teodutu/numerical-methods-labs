function [lambdas step] = QR(A, tol, max_iter)
  [n n] = size(A);
  lambdas = [];

  for step = 1:max_iter
    [Q R] = qr(A);
    A = R * Q;
    
    % verificam daca taiem sau nu
    if abs(A(2, 1)) < tol
      lambdas = [lambdas; A(1, 1)];
      A = A(2:n, 2:n);
      n--;
    endif
    
    if n == 1
      lambdas = [lambdas; A(1, 1)];
      break;
    endif
    
    if abs(A(n, n - 1)) < tol
      lambdas = [lambdas; A(n, n)];
      n--;
      A = A(1:n, 1:n);
    endif
    
    if n == 1
      lambdas = [lambdas; A(1, 1)];
      break;
    endif
    
    % daca gasim un b_j = 0 altundeva, impartim matricea
    % in 2: deasupra si dedesubtul lui b_j si rulam QR
    % pe fiecare
    % if un b_j == 0
    %  lambdas1 = QR(A(1:j, 1:j), tol, max_iter);
    %  lambdas = [lambdas; lambdas1];
    %  A = A((j+1):n, (j+1):n);
    %  n -= j;
    % endif
  endfor
endfunction