function [lambdas step] = QR_Shift(A, tol, max_iter)
  [n n] = size(A);
  lambdas = [];

  for step = 1:max_iter
    E = A((n-1):n, (n-1):n);  % coltul dreapta-jos al lui A
    eig_E = eig(E);
   
    % caut lambda' (eig_E) cel mai apropiat de A(n, n)
    if abs(eig_E(1) - A(n, n)) < abs(eig_E(2) - A(n, n))
      sigma = eig_E(1);
    else
      sigma = eig_E(2);
    endif
    
    [Q R] = qr(A - sigma * eye(n));
    A = R * Q + sigma * eye(n);  % punem sigma la loc ca sa nu stricam lambda
    
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
    % in 2: deasupra si dedesubtul lui b_j si rulam QR_Shift
    % pe fiecare
    % if un b_j == 0
    %  lambdas1 = QR_Shift(A(1:j, 1:j), tol, max_iter);
    %  lambdas = [lambdas; lambdas1];
    %  A = A((j+1):n, (j+1):n);
    %  n -= j;
    % endif
  endfor
endfunction