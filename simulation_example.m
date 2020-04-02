clear ; close all;% clc

days_to_s = 24*3600;

% simulation box 0 to L
L = 5; %[m]
% number of grid points
N_x = 150;
x = linspace(0,L,N_x)';
% simulate until
t_end = 300 *days_to_s; %[s]

% all times for which concentration profile c is computed
% use more than init and final value, since otherwise output of ode23 unexpected
t = linspace (0, t_end, 3)';



% parameters (currently all independent!)
D = 8.5 * 10^-9; % [m^2 / s]
v = 9.5 * 10^-8; % [m/s]

lambda = log(2)/ (12.32 * 365 * days_to_s); %[1/s]



% influx of concentration from surface (until t_0 surface layer has c = c_0)
% example: until t_0(1) , c is c_0(1). Then c = c_0(2) until t_0(2).
t_0 = [0.1*t_end, 0.2*t_end];
c_0 = [1, 0.2];

% initial values of concentration profile
c_init = zeros(N_x,1);






% compute concentration profiles at times t
c = transport_radio_nuclei(t, L, c_init, c_0, t_0, D, v, lambda, N_x);

% fit example
% save c profile at t_end
csvwrite('data', c(end,:)');

fit_results = fit(t, L, c_init, c_0, t_0, D, v, lambda, N_x, [D/4, v/2]);

fit_results

pause();


plot(t / days_to_s, c(:,end));
xlabel('t [days]');
ylabel('c (x=L) [a.u.]');

%
figure;

labels = {};


for i = 1:size(t,1)
  plot (x, c(i,:) );
  labels = {labels{:}, ["t = ", num2str(t(i)/days_to_s)]};
  hold on;
  %fprintf('total amount of material in box %f \n', sum(c(i,:)*L));
endfor

plot(x, csvread('data')', 'marker', 'o');

xlabel('x [m]');
ylabel('c [a.u.]');
legend(labels);