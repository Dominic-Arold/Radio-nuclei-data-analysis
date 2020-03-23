clear ; close all;% clc

days_to_s = 24*3600;



% parameters (currently all independent!)
D_L = 8.5 * 10^-9; % [m^2 / s]
v_a = 9.5 * 10^-8; % [m/s]
R_d = 1.0;
lambda = log(2)/ (12.32 * 365 * days_to_s); %[1/s]
% simulation box 0 to L
L = 1; %[m]
% number of grid points
N_x = 150;
% simulate until
t_end = 300 *days_to_s; %[s]

% all times for which concentration profile c is computed
t = linspace (0, t_end, 10)';

% influx of concentration from surface (until t_0 surface layer has c = c_0)
% example: until t_0(1) , c is c_0(1). Then c = c_0(2) until t_0(2).
t_0 = [0.1*t_end, 0.2*t_end];
c_0 = [1, 0.2];

% initial values of concentration profile
c_init = zeros(N_x,1);
%c_init(1:N_x/10) = c_0 * (1-linspace(0,1,N_x/10));

% compute concentration profiles at times t
c = transport_radio_nuclei(t, L, c_init, c_0, t_0, D_L, v_a, R_d, lambda, N_x);


plot(t / days_to_s, c(:,end));
xlabel('t [days]');
ylabel('c (x=L) [a.u.]');

%pause();
figure;

labels = {};
x = linspace(0,L,N_x)';

for i = 1:size(t,1)
  plot (x, c(i,:) );
  labels = {labels{:}, ["t = ", num2str(t(i)/days_to_s)]};
  hold on;
  %fprintf('total amount of material in box %f \n', sum(c(i,:)*L));
endfor
xlabel('x [m]');
ylabel('c [a.u.]');
legend(labels);