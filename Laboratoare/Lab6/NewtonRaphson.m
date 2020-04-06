function [x, i] = NewtonRaphson(f, df, x0, tol, max_iter)
  for i = 1:max_iter
    x = x0 - f(x0) / df(x0);
    
    if abs(f(x)) < eps
      break;
    endif

    if abs(x - x0) < tol
      break;
    endif

    x0 = x;
  endfor
endfunction