clear ; close all;% clc
format short e;

% constant for time cenversion
days_to_s = 24*3600;

% simulation box 0 to L in meter
L = 4; %[m]
% number of grid points
N_x = 150;
% simulate until
t_end = 400 *days_to_s; %[s]
% number of time points for which to return simulations results
N_t = 10;

% all times in vector t for which concentration profile c will be computed
% (use more than init and final value, since otherwise output of ode23 solver not handled correctly)
t = linspace(0, t_end, N_t)';
x = linspace(0,L,N_x)';


% polinomial coefficients of spacially varying parameters D, v 
% In this example D is D0 at surface x=0 and linearly increases to 1.75*D0 at x=L.
% same for v or v is given from D via the dispersion length (v = D/L_D).
D0 = 8.5 * 10^-9; % [m^2 / s]
D = [D0, 0.75*D0 / L];
%v0 =  9.5 * 10^-8; % [m/s]
%v = [v0, 0.5*v0 / L];
v = D / 0.1;
% decay constant
lambda = log(2)/ (12.32 * 365 * days_to_s); %lambda (Tritium) [1/s]


% ===============================================================
% Example 1: Computed concentration profiles.

% influx of concentration from surface (until t_0 surface layer has c = c_0)
% until t_0(1) , c is c_0(1). Then c = c_0(2) until t_0(2), etc.
t_0 = [0.1*t_end, 0.2*t_end];
c_0 = [1, 1];

% initial values of concentration profile at t = 0
c_init = zeros(N_x,1);


% compute concentration profiles at times t
c = transport_radio_nuclei(t, L, c_init, c_0, t_0, D, v, lambda, N_x);


plot(t / days_to_s, c(:,end));
xlabel('t [days]');
ylabel('c (x=L) [a.u.]');
title('Example 1: Concentration at end of simulation box.');

figure;

labels = {};
for i = 1:N_t
  plot (x, c(i,:) );
  labels = {labels{:}, ["t = ", num2str(t(i)/days_to_s)]};
  hold on;
  %fprintf('total amount of material in box %f \n', sum(c(i,:)*L/N_x));
endfor

xlabel('x [m]');
ylabel('c [a.u.]');
legend(labels);
title('Example 1: Computed concentration profiles with linearly varying D,v.');




% ===============================================================
% Example 2: Fit to artificial data set with D,v once assumed constant and once linearly varying in space.
fprintf('Starting fit example.\n');
pkg load optim;

% save c profile at t_end with additinal noise to imitate experimental data for fit example
noise_range = 0.02;
%point_dist = int32(N_x/ 33);
N_data = 25;
prefix = 'data_t_';

% add noise and save
for i = 2:N_t
  fname = strcat(prefix, num2str(t(i)), '.csv');
  c_data = c(i,:)';
  noise = noise_range*(2*rand(N_x,1)-1);
  c_data += c_data .* noise;
  dataset = [x,c_data];
  csvwrite(fname, dataset(1:int32(N_x/ N_data):end,:));
endfor

% initial parameter guess of D,v for fit
D_init = [0.9*D(1),0];
v_init = [0.8*v(1),0];
% time t indices from which data is used in the fit
i_fit = [5,7,10];
% corresponding file names in cell array
fnames = {};
for i = 1:max(size(i_fit))
  fnames{end+1} = strcat(prefix, num2str(t(i_fit(i))), '.csv');
endfor


% perform fit to noisy data 
[p_const,f] = fit(fnames, t(i_fit), L, c_init, c_0, t_0, D_init(1), v_init(1), lambda, N_x);
fprintf('Got coefficients for const. D,v:\n');
disp(p_const);

[p,f] = fit(fnames, t(i_fit), L, c_init, c_0, t_0, D_init, v_init, lambda, N_x);
fprintf('Got coefficients for linearly varying D,v:\n');
disp(p);

fprintf('Originally used:\n');
disp(D');
disp(v');


% plot saved datasets, and optimized c profiles for constant and linearly varying...
% parameters [D,v], at all times used in the optimization.
order_D = max(size(D));
c_opt = transport_radio_nuclei(t(i_fit), L, c_init, c_0, t_0, p(1:order_D), ...
                                p(order_D+1:end), lambda, N_x)(2:end,:);
c_opt_const = transport_radio_nuclei(t(i_fit), L, c_init, c_0, t_0, p_const(1), ...
                                p_const(2), lambda, N_x)(2:end,:);

figure;                           
labels = {};
for i = 1:size(i_fit,2)
  %plot(x, c(i_fit(i),:), 'linestyle', '--');
  %hold on;
  dataset = csvread(fnames{i});
  plot(dataset(:,1), dataset(:,2) , 'marker', 'o', 'linestyle', 'none');
  labels{end+1} = strcat("data  t = ", num2str(t(i_fit(i))/days_to_s));
  hold on;
  
  plot(x, c_opt_const(i,:), 'linestyle', '--');
  labels{end+1} = strcat("fit (const)  t = ", num2str(t(i_fit(i))/days_to_s));
  hold on;
  
  plot(x, c_opt(i,:));
  labels{end+1} = strcat("fit (1.order) t = ", num2str(t(i_fit(i))/days_to_s));
  hold on;
endfor
xlabel('x [m]');
ylabel('c [a.u.]');
legend(labels);
title('Example 2: Fit to artificial data set with D,v once assumed constant and once linearly varying in space.');

% ===============================================================