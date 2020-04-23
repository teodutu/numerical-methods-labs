function [lambda, v, step] = MPD(A, v, tol, max_iter)
  lambda = Inf;
  
  % v != 0
  
  for step = 1:max_iter
    z = A * v;
    % daca v = 0 => z = 0 => v = 0 / 0
    v = z / norm(z);  % ca sa scapam de a1
    
    lambda_new = v' * A * v;
    
    if norm(lambda_new - lambda) < tol  % lambda e complex, deci folosim norma
      lambda = lambda_new;
      break;
    endif
    
    lambda = lambda_new;
  endfor
endfunction