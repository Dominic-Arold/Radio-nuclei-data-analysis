function xdot = f( t, x, cur_coef, dx, R_d, t_0, c_0, v)

xdot = cur_coef * x;
xdot(1) += J_1_in(t, t_0, c_0, v) / (dx*R_d);

endfunction