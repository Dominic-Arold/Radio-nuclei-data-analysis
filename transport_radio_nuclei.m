function result = transport_radio_nuclei(t, L, c_init, c_0, t_0, ...
    D, v, lambda, N_x, return_all_result = 1, rescaling = 1)
%{
Use the transport equation for radio nuclei to compute the concentration profile c...
at different times.
  
Returns: 2 dim. array where every row is the concentration profile at a given time point in t.        

Input arguments:
t                   : (1 dim. array) List of times for which concentration profile c ...
                      will be computed. First entry needs to be 0, otherwise added.
c_0                 : (1 dim. array) List of inflowing concentration values ...
                      at surface until time in t_0.
t_0                 : (1 dim. array like c_0) List of times until which inflowing ...
                       concentration c_0 persists.
D, v                : (1 dim. array) Lists of Taylor coefficients of spatially varying ...
                       diffusion and drift velocity in transport equation.
lambda              : (float) Decay constant of radio nuclide.
L                   : (float) length of the simulation box.
c_init              : (N_x x 1 array) Initial distribution of concentration ...
                       profile at t = 0.
return_all_results  : (bool) Return all results or just the ...
                      concentration profile c of the last time in t.
rescaling           : (bool) Switching on/off rescaling for transport equation.
%}


R_d = 1.0;

order_D = max(size(D));
order_v = max(size(v));
order_Dv = min(order_D, order_v);

if(order_D != order_v)
  warning('Polynomial order of diffusion and drift coefficient, %i and %i , are assumed to be equal. Continue with order %i for both.\n',...
  order_D, order_v, order_Dv);
endif

if(size(t,1) != 1)
  t = t';
endif
if(t(1) != 0)
  t = [0,t];
endif

% ===============================================================
% rescaling
if(rescaling)
  X = L; % m
  T = t(end); % s
  C = 1; % 
  
  t /= T;
  L /= X;
  c_init /= C;
  c_0 /= C;
  t_0 /= T;
  for i = 1:order_Dv
    D(i) /= (X^2/T) * 1/X^(i-1);
    v(i) /= (X/T) * 1/X^(i-1);
  endfor
  
  lambda /= (1.0/T);
endif
% ===============================================================

dx = L/N_x;
x = linspace(0,L,N_x)';

dummy = zeros(2,N_x);
for i = 1:N_x
  for j = 1:order_Dv
    dummy(1,i) += D(j) * x(i)^(j-1);
    dummy(2,i) += v(j) * x(i)^(j-1);
  endfor
endfor
D = dummy(1,:);
v = dummy(2,:);


% ===============================================================
% define current matrix
% J_1_in not included
cur_coef = zeros(N_x,N_x);

% J_in und out for 2 to N-1
for i = 2:N_x-1
  cur_coef(i, i-1) = (D(i-1)+D(i))/2.0/dx + (v(i-1)+v(i))/2.0/2.0;
  cur_coef(i, i) = -(D(i-1)+2*D(i)+D(i+1))/2.0/dx + (v(i-1)-v(i+1))/2.0/2.0;
  cur_coef(i, i+1) = (D(i)+D(i+1))/2.0/dx - (v(i)+v(i+1))/2.0/2.0;
endfor

% J_1_out
cur_coef(1,1) -=  (D(1)+D(2))/2.0/dx + (v(1)+v(2))/2.0/2.0;
cur_coef(1,2) -= -(D(1)+D(2))/2.0/dx + (v(1)+v(2))/2.0/2.0;

% J_N_in (is J_(N-1)_out)
cur_coef(N_x,N_x-1) += (D(N_x-1)+D(N_x))/2.0/dx + (v(N_x-1)+v(N_x))/2.0/2.0;
cur_coef(N_x,N_x) +=  -(D(N_x-1)+D(N_x))/2.0/dx + (v(N_x-1)+v(N_x))/2.0/2.0;

% J_N_out
cur_coef(N_x,N_x-1) -=  -(v(N_x-1)+v(N_x))/2.0 / 2.0;
cur_coef(N_x,N_x) -=  (v(N_x-1)+v(N_x))/2.0 * 1.5;


cur_coef = cur_coef / (R_d*dx);

for i = 1:N_x
  cur_coef(i,i) -= lambda;
endfor


% ===============================================================
% time iteration
opt = odeset ("RelTol", 10^(-6), "AbsTol", 10^(-6));
cdot = @(t,c)(f(t, c, cur_coef, dx, R_d, t_0, c_0, v(1)));

start = time();
% cdot is integrated with the command
[t,c] = ode23(cdot, t, c_init, opt);
%fprintf('Execution time: %f s\n', time()-start);




if (return_all_result)
  result = c;
else
  result = c(end,:);
endif

endfunction