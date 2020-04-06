function [x, i] = Bisection(f, a, b, tol, max_iter)
  if f(a) * f(b) > 0
    x = NaN;
    i = -1;
  endif
  
  x_prev = inf;
  
  for i = 1:max_iter
    x = (a + b) / 2

    if abs(x - x_prev) < tol
      break;
    endif

    p = f(a) * f(x);

    if p < 0
      b = x;
    elseif p > 0
      a = x;
    else
      break;
    endif
    
    x_prev = x;
  endfor
endfunction