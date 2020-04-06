function [x, i] = SecantMethod(f, x0, x1, tol, max_iter)
  for i = 1:max_iter
    f_x1 = f(x1);  % f(x^(p))
    f_x0 = f(x0);  % f(x^(p - 1))

    x = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
    
    if abs(f(x)) < eps
      break;
    endif

    if abs(x - x1) < tol
      break;
    endif

    x0 = x1;
    x1 = x;
  endfor
endfunction