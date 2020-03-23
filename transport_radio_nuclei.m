function c = transport_radio_nuclei(t, L, c_init, c_0, t_0, ...
  D_L, v_a, R_d, lambda, N_x, step_size_scale = 5.0)



% ===============================================================
% rescaling
X = L; % m
T = 100 * 24*3600; % s
C = 1; % 

t /= T;
L /= X;
c_init /= C;
c_0 /= C;
t_0 /= T;
D_L /= (X^2/T);
v_a /= (X/T);
R_d /= 1.0;
lambda /= (1.0/T);

% ===============================================================


dx = L/N_x;
%dt = dx / v_a / step_size_scale;
%N_t = uint64(t(end) / dt);



% ===============================================================
% define current matrix
% J_1_in not included
cur_coef = zeros(N_x,N_x);

% J_in und out for 2 to N-1
for i = 2:N_x-1
  cur_coef(i, i-1) = (D_L/dx + v_a/2.0);
  cur_coef(i, i) = (-2*D_L/dx);
  cur_coef(i, i+1) = (D_L/dx - v_a/2.0);
endfor

% J_1_out
cur_coef(1,1) -= (D_L/dx + v_a/2.0);
cur_coef(1,2) -= (-D_L/dx + v_a/2.0);

% J_N_in (is J_(N-1)_out)
cur_coef(N_x,N_x-1) += (D_L/dx + v_a/2.0);
cur_coef(N_x,N_x) += (-D_L/dx + v_a/2.0);

% J_N_out
cur_coef(N_x,N_x-1) -=  -v_a / 2.0;
cur_coef(N_x,N_x) -=  v_a * 1.5;


cur_coef = cur_coef / (R_d*dx);

for i = 1:N_x
  cur_coef(i,i) -= lambda;
endfor


% ===============================================================

cdot = @(c,t)(f(c, t, cur_coef, dx, R_d, t_0, c_0, v_a));
% cdot is integrated with the command
%[c, istate, msg] = lsode(cdot, c_init, t);
[t,c] = ode45(cdot, t, c_init);










%{
r_count = 1;

c = c_init;
for i = 0:N_t
  dc = cdot(c, i*dt) * dt;  
  
  DC(i+1) = sum(dc);
  if any(0:9 == i) %any(uint64(t/dt) == i) 
    results(r_count,:) = c';
    r_count += 1;  
  endif

  c += dc;
endfor
%}

%fprintf('Change in c: %f \n', sum(DC));
%plot(linspace(1,N_t,size(DC,2)), DC);



endfunction