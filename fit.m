function result = fit(t, L, c_init, c_0, t_0, ...
      D, v, lambda, N_x, p_init)


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
D /= (X^2/T);
v /= (X/T);
%R_d /= 1.0;
lambda /= (1.0/T);

p_init(1) /= (X^2/T);
p_init(2) /= (X/T);

% ===============================================================

  

fitfunc = @(dummy, p)(transport_radio_nuclei(t, L, c_init, c_0, t_0, p(1), p(2),...
                lambda, N_x, 1.0, 0, 0));
                
[f, p, cvg, iter, corp, covp, covr, stdresid, Z, r2] = ...
leasqr([], csvread('data')', p_init, fitfunc);

D = p(1) * (X^2/T);
v = p(2) * (X/T);


result = [D,v];
endfunction
