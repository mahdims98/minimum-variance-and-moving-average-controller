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

%% minimum variance_not_adaptive.
rng(50);

syms q;
main_folder_name = "1-4-1-MV";
main_folder = 'images/q1-final/' + main_folder_name + "/";
sub_name = "non_adaptive_MV_";
% main_title = main_folder + sub_name;


% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;

% noise and its degree
% noise_poly = [1,0.4,0.2].';
noise_poly = [1,0,0].';

C = noise_poly.';

deg_noise = length(noise_poly);
noise_variance = 0.005; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
% v = [zeros([1,ceil(num_samples/2)]), 10*ones([1,ceil(num_samples/2)])];
v = zeros([1,num_samples]);


% must be changed for every problem
len_desA = 3;
len_desB = 2;
deg_desA = len_desA - 1;
deg_desB = len_desB - 1;
% A_estimated = [1,0,0]; % initial conditions
% B_estimated = [0.01, 0.035];
% assert(length(A_estimated)==deg_desA && length(B_estimated)==deg_desB, "initial condistions are not true")

% *********************
% MV parameter
% *********************

%  [Bplus, Bminus] = factor_polynomial(B);
% d = deg_desA - (length(Bplus)-1); % deg_desA - deg_Bplus

d0 = deg_desA - deg_desB;
% manually not cancelling the zeros
d = 2;
Bplus = 1;

assert(d>=d0, "d < d0 !")

A_dio = A;
B_dio = B; 
D_dio = conv([1 , zeros([1,d-1])], C); % RHS of the diophantin eqn.
D_dio = conv(D_dio, Bplus);
[alpha, beta] = solve_diophantin_general(A_dio, B_dio, D_dio, 0); % q^(d-1)*C*Bplus = A*alpha+B*beta // (F=Alpha, G=Beta)
R = alpha; 
S = beta; 

% initial conditions for y
skip_instances = max(len_desA, len_desB);
total_parameters = len_desA - 1 + len_desB;
y = [];
y(1:skip_instances) = 0.1;
u(1:skip_instances) = 0;

theta_real = [A(2:end).'; B.'];

% for plotting
theta_hat_toPlot = zeros([num_samples, total_parameters]);

for i = skip_instances:num_samples
    
    phi_t = [-y(i-1:-1:i-(len_desA - 1)), u(i-1:-1:i-len_desB)].';  

    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly;
    y(i) = phi_t.' * theta_real + noise_t + B * v(i:-1:i-(length(B)-1)).';
    
    u(i) = (S * [-y(i:-1:i-(length(S)-1))].' - R(2:end) * [u(i-1:-1:i-(length(R)-1))].')/R(1);

end

%% metrics
varY = var(y);
varU = var(u);
varabsY = var(abs(y));
varabsU = var(abs(u));
meanY = mean(y);
meanU = mean(u);

acclosY = zeros([1, length(y)]);
acclosY(1) = y(1)^2;

for i=2:length(y)
    acclosY(i) = y(i)^2 + acclosY(i-1);
end

acclosU = zeros([1, length(u)]);
acclosU(1) = u(1)^2;

for i=2:length(u)
    acclosU(i) = u(i)^2 + acclosU(i-1);
end

%% plotters

if ~exist(main_folder, 'dir')
   mkdir(main_folder)
end

% writing parameters
metrics_matrxi = [varY, varU, varabsY, varabsU, meanY, meanU]; 
metrics_table = array2table(metrics_matrxi);
metrics_table.Properties.VariableNames(1:end) = {'varY','varU','varabsY', 'varabsU', 'meanY', 'meanU'};
writetable(metrics_table,main_folder + sub_name + '.csv')

figure()
subplot(2,1,1);
plot(t,y, 'DisplayName','Real')
hold on;

xlabel("sample number");
title("output");
legend('Location','best');
subplot(2,1,2);
plot(t, u, 'DisplayName','control input')
xlabel("sample number");
title("control input");
saveas(gcf, main_folder + sub_name + "y-u" + '.jpeg')
close all

% metrics
f6 = figure();
subplot(2,1,1);
plot(t,acclosY,'DisplayName',' output accumulated loss')
% hold on; 
% plot(t,acclosY_indirect, 'DisplayName','indirect')
legend();
title("y accumulated loss");
subplot(2,1,2);
plot(t,acclosU,'DisplayName','control accumulated loss')
title("u accumulated loss");
saveas(gcf, main_folder + sub_name + "accu-loss" + '.jpeg')
close all

%%
((A(3)-C(3))*B(1)*B(2) + (C(2)-A(2))*B(2)^2)/(B(2)^2+A(2)*B(1)*B(2)+A(3)*B(1)^2)
(-B(2) * (A(2)^2-A(3)-C(2)*A(2)+C(3))+B(1)*(C(2)*A(3)-A(2)*A(3)))/(B(2)^2+A(2)*B(1)*B(2)+A(3)*B(1)^2)
(-B(2) * (A(2)*A(3)-C(2)*A(3)) + B(1)*(A(3)*C(3)-A(3)^2))/(B(2)^2+A(2)*B(1)*B(2)+A(3)*B(1)^2)