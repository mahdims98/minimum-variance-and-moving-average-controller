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

%% indirect - overparameter
rng(50);

main_folder_name = "1-8-1-minVar-underparameter";
main_folder = 'images/q1-final/' + main_folder_name + "/";
sub_name = "overparameter";
main_title = sub_name;
close_all = true;

% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;

% noise and its degree
% noise_poly = [1,0.4,0.2].';
noise_poly = [1,0.4,0.2].';
deg_noise = length(noise_poly);
noise_variance = 0.005; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
% v = [zeros([1,ceil(num_samples/2)]), 10*ones([1,ceil(num_samples/2)])];
v = zeros([1,num_samples]);



% must be changed for every problem
% initial conditions
len_desA = 2;
len_desB = 1;
deg_desA = len_desA - 1;
deg_desB = len_desB;
A_estimated = [1,0]; % initial conditions
B_estimated = [0.01];
% RLS initial parameters (C)
theta_epsilon_zero = [0.1];
C_estimated = [1, theta_epsilon_zero.']; % C must be monic 

assert(length(A_estimated)==len_desA && length(B_estimated)==len_desB, "initial condistions are not true")

d0 = deg_desA - deg_desB;

% initial R S T and Ao
R_calculated = [1]; 
S_calculated = [1];
 

% R S which were calculated in the last part
R_real = [0.4200, -0.1890];
S_real = [1.0700, 0.0950];



% initial conditions for y
skip_instances = max(length(A), length(B));
total_parameters = len_desA - 1 + len_desB;
y = [];
y(1:skip_instances) = 0.1;
u(1:skip_instances) = 0;

% actual parameters
theta_real = [A(2:end).'; B.'];
theta_real_toplot = [A(2:end).'; B.'; noise_poly(2:end)];


% for plotting
% u_toPlot = zeros([num_samples, length(uc)]);
R_calculated_toPlot = zeros([num_samples, length(R_calculated)]);
S_calculated_toPlot = zeros([num_samples, length(S_calculated)]);


theta_hat_toPlot = zeros([num_samples, len_desA - 1 + len_desB + 1]);

els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);
for i = skip_instances:num_samples
    phi_t_real = [-y(i-1:-1:i-(length(A) - 1)), u(i-1:-1:i-length(B))].';
    phi_t = [-y(i-1:-1:i-(len_desA - 1)), u(i-1:-1:i-len_desB)].';  
    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly;
    y(i) = phi_t_real.' * theta_real + noise_t + B * v(i:-1:i-(length(B)-1)).';
    
    u(i) = (S_calculated * [-y(i:-1:i-(length(S_calculated)-1))].' - R_calculated(2:end) * [u(i-1:-1:i-(length(R_calculated)-1))].') /R_calculated(1);
    theta_hat_new = els_solver.update_ELS(y(i), phi_t);

    A_estimated = [1, theta_hat_new(1:(len_desA - 1)).'];
    B_estimated = theta_hat_new(len_desA:total_parameters).';
    C_not_modified = [1, theta_hat_new(total_parameters+1:end).'];
    C_estimated = modify_polynomial(C_not_modified); % removing zeros outside unit circle
    
    % controller parameters
    A_dio = A_estimated;
    B_dio = B_estimated; 
    D_dio = conv([1 , zeros([1,d0-1])], C_estimated); % RHS of the diophantin eqn.
    D_dio = conv(D_dio, B_estimated);
    [alpha, beta] = solve_diophantin_general(A_dio, B_dio, D_dio, 0); % q^(d0-1)*C*Bplus = A*alpha+B*beta // (F=Alpha, G=Beta)
    R_calculated = alpha; 
    S_calculated = beta; 

    if R_calculated(1) == 0
        R_calculated = R_calculated(2:end);
    end
    if S_calculated(1) == 0
        S_calculated = S_calculated(2:end);
    end

    theta_hat_toPlot(i, :) = theta_hat_new(1:end).';
    R_calculated_toPlot(i, :) = R_calculated;
    S_calculated_toPlot(i, :) = S_calculated;

%     disp(R_solved)
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

% acclosY_indirect = acclosY
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
xlabel("sample number");
title("output");
legend('Location','best');
subplot(2,1,2);
plot(t, u, 'DisplayName','control input')
xlabel("sample number");
title("control input");
saveas(gcf, main_folder + sub_name + "y-u" + '.jpeg')


% metrics
f6 = figure();
subplot(2,1,1);
plot(t,acclosY,'DisplayName',' output accumulated loss')
title("y accumulated loss");
subplot(2,1,2);
plot(t,acclosU,'DisplayName','control accumulated loss')
title("u accumulated loss");
saveas(gcf, main_folder + sub_name + "accu-loss" + '.jpeg')


% R 
total_plot_rows = ceil(length(R_calculated)/2) + ceil(length(S_calculated)/2)-1;
f2 = figure();
for i = 1:(length(R_calculated)) % first parameter of R is 1
    subplot(total_plot_rows,2,i);
    title_text = "R_%d";
    plot(1:num_samples, R_calculated_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")

end
end_of_r_plot = i;

% S
for i = 1:(length(S_calculated)) % first parameter of R is 1
    subplot(total_plot_rows,2,i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_calculated_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,main_folder + sub_name + "-S-R" + '.jpeg')
% close all

% Theta
f5 = figure();
f5.Position = [0,0,400,900];
for i = 1:length(theta_hat_toPlot(1,:))
    title_text = "θ_%d";
    subplot(ceil(length(theta_hat_toPlot(1,:))/2),2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'DisplayName','Predicted')
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * theta_real_toplot(i) , 'DisplayName','Real')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,main_folder + sub_name + "-theta" + '.jpeg')
if close_all
    close all
end