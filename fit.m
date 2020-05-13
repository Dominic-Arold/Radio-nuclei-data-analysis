function [params,f] = fit(fnames, t, L, c_init, c_0, t_0, ...
      D, v, lambda, N_x, use_L_D = false)

%{
Determine spatially varying diffusion and drift velocity coefficients by performing...
a least square optimization of the transport equation to given concentration...
profiles at several times.

Returns:
params : (1 dim. array) Taylor coefficients of spatially varying diffusion and...
          drift velocity determined from the optimization. First half of ...
          parameters corresponds to D and the second to v, meaning...
          params = [D_opt, v_opt].
f      : (2 dim. array) Every row is the found optimized concentration profile...
          at a given time point in t corresponding to [D,v] = params. Only points ...
          corresponding to measurement positions at each time are returned.      

Input arguments:
fnames    : (1 dim. cell array) List of file name strings. Each file has to contain ...
             two columns with first the spatial positions of measurement and ...
             second the measured concentration profile used in the optimization. ...
             All files need to hve the same number of rows.
t         : (1 dim. array like fnames) List of times corresponding to file data ...
             given by fnames.
L         : (float) length of the simulation box used for optimization. Needs ...
             to contain all positions listed in the input files fnames.
c_init    : (1 dim. array) Initial concentration profile at t = 0.
c_0       : (1 dim. array) List of inflowing concentration values ...
             at surface until time in t_0.
t_0       : (1 dim. array like c_0) List of times until which inflowing ...
             concentration c_0 persists.
D, v      : (1 dim. array) Lists of initial guesses for Taylor coefficients of...
             spatially varying diffusion and drift velocity coefficients. ...
             When use_L_D is true, D and v must have same length, i.e. number ...
             of orders considered. Otherwise, v must be float.
lambda    : (float) Decay constant of radio nuclide.
N_x       : (integer) number of spacial grid points. Should be chosen high ...
             enough to have a higher grid point density than measurement point ...
             density.
use_L_D   : (bool) When false, all Taylor coefficients D,v are independent fit ...
             parameters. When true, the coefficients D are fit paramters and ...
             the drift velocity coefficients are determined as v = D / L_D ...
             with L_D only one additional fit parameter.
%}
      
      

order_D = max(size(D));
order_v = max(size(v));

if(use_L_D)
  if(order_v != 1)
    error('When "use_L_D", the input variable v must be scalar initial value.');
  endif
elseif(order_D != order_v)
  error('Polynomial order of diffusion and drift coefficient, %i and %i , are assumed to be equal.\n',...
  order_D, order_v);
endif


% ===============================================================
% rescaling
  X = L; % m
  T = t(end); % s
  C = 1; % 
  
  t /= T;
  L /= X;
  c_init /= C;
  c_0 /= C;
  t_0 /= T;
  for i = 1:order_D
    D(i) /= (X^2/T) * 1/X^(i-1);
  endfor
  for i = 1:order_v
    v(i) /= (X/T) * 1/X^(i-1);
  endfor
  lambda /= (1.0/T);
  
  
% ===============================================================
% load file data
for i = 1:max(size(t))
  filedata = csvread(fnames{i});
  x_data(i,:) = filedata(:,1)' / X;
  c_data(i,:) = filedata(:,2)' / C;
endfor


x = linspace(0,L,N_x);

if(size(t,1) != 1)
  t = t';
endif
if(t(1) != 0)
  t = [0,t];
endif

% ===============================================================
% fit function
core_func = @(D,v)(filter_to_exp_data(transport_radio_nuclei(
                  t, L, c_init, c_0, t_0, D, v, lambda, N_x, 1, 0)(2:end,:),
                  x, x_data));
if(use_L_D)
  p_init = cat(2,D, D(1)/v);
  fitfunc = @(dummy, p)(core_func(p(1:order_D), p(1:order_D)/p(order_D+1)));
else
  p_init = cat(2,D,v);
  fitfunc = @(dummy, p)(core_func(p(1:order_D), p(order_D+1:2*order_D)));
endif


% ===============================================================
% perform optimization
[f, p, cvg, iter, corp, covp, covr, stdresid, Z, r2] = ...
leasqr([], c_data, p_init, fitfunc);

fprintf('Fit finished after %i iterations.\n', iter);


% ===============================================================
% format output
if(use_L_D)
  for i = 1:order_D
    p(i) *= (X^2/T) * 1/X^(i-1);
  endfor
  p(order_D+1) *= X;
  params = cat(1,p(1:order_D), p(1:order_D)/p(order_D+1));
else
  for i = 1:order_D
    p(i) *= (X^2/T) * 1/X^(i-1);
    p(order_D+i) *= (X/T) * 1/X^(i-1);
  endfor
  params = p;
endif


endfunction
