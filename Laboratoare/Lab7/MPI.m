function [lambda, v, step] = MPI(A, v, u, tol, max_iter)
  lambda = Inf;
  
  [n n] = size(A);
  
  for step = 1:max_iter
    z = (A - u * eye(n)) \ v;
    v = z / norm(z);
    
    lambda_new = v' * A * v;
    
    if norm(lambda_new - lambda) < tol  % lambda e complex, deci folosim norma
      lambda = lambda_new;
      break;
    endif
    
    lambda = lambda_new;
  endfor
endfunction