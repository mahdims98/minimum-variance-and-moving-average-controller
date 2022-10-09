%% discretization of the model
clear;

Ts = 0.2;
sysD = zpk(0.45,[0.25,0.42],0.42,Ts);
[numD, denD] = tfdata(sysD, 'v');


% bodeplot(sysD);

B = numD;
A = denD;
% step(sysD)
%%
if B(1) == 0
    B = B(2:end);
end
%% desired system

overshoot = 20;
rise_time = 0.3;

zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 1.8/(rise_time); 

% sys_2d_desired  = tf([wn^2], [1, 2*zeta*wn, wn^2]);
% damp(sys_desired)
z1 = -15;

k2 = -1/z1;

G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1],[],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

BmPrime_main = num_desD; 
Am = den_desD;

%%
if B(1) == 0
    B = B(2:end);
end
if BmPrime_main(1) == 0
    BmPrime_main = BmPrime_main(2:end);
end
%% minimum variance_not_adaptive
syms q;
main_folder_name = "1-4-1-min-var";
main_folder = 'images/q1/' + main_folder_name + "/";
sub_name = "non_adaptivemin_var";
% main_title = main_folder + sub_name;


% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
input_noise_variance = 0.001;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);
uc(1) = -1;
uc = uc + input_noise;

% noise and its degree
% noise_poly = [1,0.4,0.2].';
noise_poly = [1,0,0].';
C = noise_poly;
deg_noise = length(noise_poly);
noise_variance = 0.0; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
% v = [zeros([1,ceil(num_samples/2)]), 10*ones([1,ceil(num_samples/2)])];
v = zeros([1,num_samples]);

% must be changed for every problem
deg_desA = 3;
deg_desB = 2;
A_estimated = [1,0,0]; % initial conditions
B_estimated = [0.01, 0.035];

assert(length(A_estimated)==deg_desA && length(B_estimated)==deg_desB, "initial condistions are not true")

% solving diophantin
d0 = degA - degB;
[alpha, beta] = solve_diophantin_general(A, B, D, delay)

% R and S



% initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;

theta_real = [A(2:end).'; B.'];

% for plotting
theta_hat_toPlot = zeros([num_samples, total_parameters]);

% RLS initial parameters
theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);

for i = skip_instances:length(uc)
    
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';  

    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly;
    y(i) = phi_t.' * theta_real + noise_t + B * v(i:-1:i-(length(B)-1)).';
    
    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_solved * [-y(i:-1:i-(length(S_solved)-1))].' - R_solved(2:end) * [u(i-1:-1:i-(length(R_solved)-1))].';

    theta_hat_new = els_solver.update_ELS(y(i), phi_t);

    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    Bplus = [B_estimated/B_estimated(1)];
    Bminus = B_estimated(1);
    BmPrime = BmPrime_main/B_estimated(1);
    Ao = [1];

    [R_solved,S_solved,T] = solve_diophantin(Am, [1,A_estimated], Ao, Bminus, Bplus, BmPrime);

    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_solved_toPlot(i, :) = R_solved;
    S_solved_toPlot(i, :) = S_solved;
    T_toPlot(i, :) = T;
%     disp(R_solved)
end
