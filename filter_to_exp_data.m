function result = filter_to_exp_data(c_num, x_num, x_data)

N_t = size(c_num,1);
N_x = max(size(x_num)); 
result = [];

for h = 1:N_t
  count = 1;
  c_filtered = [];
  for i = 1:N_x
    if(x_num(i) >= x_data(h,count))
      c_filtered = [c_filtered, c_num(h,i)];
      if(count == size(x_data,2))
        break;
      else
        count += 1;
      endif
    endif
  endfor
  result(h,:) = c_filtered;
endfor
  
endfunction
