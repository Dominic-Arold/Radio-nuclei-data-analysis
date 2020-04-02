function x = J_1_in(t, t_0, c_0, v)
  
  if t < t_0(end)
    x = v * c_0(find(t_0 > t)(1));
  else
    x = 0;
  endif
  
endfunction
