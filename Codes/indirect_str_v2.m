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


%% indirect with zero cancellation 
main_folder_name = "1-1-indirect-str";
main_folder = 'images/q1-final/' + main_folder_name + "/";
sub_name = "colored-noise";
main_title = sub_name;

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
noise_poly = [1,0.4,0.2].';
% noise_poly = [1,0,0].';
deg_noise = length(noise_poly);
noise_variance = 0.005; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
% v = [zeros([1,ceil(num_samples/2)]), 10*ones([1,ceil(num_samples/2)])];
v = zeros([1,num_samples]);



% must be changed for every problem
% initial conditions
len_desA = 3;
len_desB = 2;
deg_desA = len_desA - 1;
deg_desB = len_desB;
A_estimated = [1,0,0]; % initial conditions
B_estimated = [0.01, 0.035];

% RLS initial parameters (C)
theta_epsilon_zero = [0.1,0.1].';
C_estimated = [1, theta_epsilon_zero.']; % C must be monic 

assert(length(A_estimated)==len_desA && length(B_estimated)==len_desB, "initial condistions are not true")

% initial B
Bplus = [B_estimated/B_estimated(1)];
Bminus = B_estimated(1);
BmPrime = BmPrime_main/B_estimated(1);


% initial R S T and Ao
R_solved = [1,-0.1]; 
S_solved = [1,0.1];
T = [1.1,0];
Ao = [1];
 

% initial conditions for y
skip_instances = max(len_desA, len_desB);
total_parameters = len_desA - 1 + len_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;

% actual parameters
theta_real = [A(2:end).'; B.'];



% for plotting
R_solved_toPlot = zeros([num_samples, length(R_solved)]);
S_solved_toPlot = zeros([num_samples, length(S_solved)]);
T_toPlot = zeros([num_samples, length(T)]);


theta_hat_toPlot = zeros([num_samples, total_parameters]);

% RLS initial parameters
theta_epsilon_zero = [0.5, 0.5].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);
for i = skip_instances:num_samples

    phi_t = [-y(i-1:-1:i-(len_desA - 1)), u(i-1:-1:i-len_desB)].'; 

    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly;
    y(i) = phi_t.' * theta_real + noise_t + B * v(i:-1:i-(length(B)-1)).';
    
    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_solved * [-y(i:-1:i-(length(S_solved)-1))].' - R_solved(2:end) * [u(i-1:-1:i-(length(R_solved)-1))].';

    theta_hat_new = els_solver.update_ELS(y(i), phi_t);

    A_estimated = theta_hat_new(1:(len_desA - 1)).';
    B_estimated = theta_hat_new(len_desA:total_parameters).';
    
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
hold on;
plot(t, uc, "--r", 'DisplayName','desired')
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
title("y accumulated loss");
subplot(2,1,2);
plot(t,acclosU,'DisplayName','control accumulated loss')
title("u accumulated loss");
saveas(gcf, main_folder + sub_name + "accu-loss" + '.jpeg')
close all

% % R 
% total_plot_rows = ceil(length(R_calculated)/2) + ceil(length(S_calculated)/2);
% f2 = figure();
% f2.Position = [0 0 400 900];
% for i = 1:(length(R_calculated)) % first parameter of R is 1
%     subplot(total_plot_rows,2,i);
%     title_text = "R_%d";
%     plot(1:num_samples, R_calculated_toPlot(:,i), 'DisplayName','Predicted') 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
%     title(sprintf(title_text, i));
%     legend('Location','best');
%     xlabel("sample number")
% 
% end
% end_of_r_plot = i;
% 
% % S
% for i = 1:(length(S_calculated)) % first parameter of R is 1
%     subplot(total_plot_rows,2,i + end_of_r_plot);
%     title_text = "S_%d";
%     plot(1:num_samples, S_calculated_toPlot(:,i), 'DisplayName','Predicted') 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
%     title(sprintf(title_text, i));
%     legend('Location','best');
%     xlabel("sample number")
% end
% saveas(gcf,main_folder + sub_name + "-S-R" + '.jpeg')
% % close all

% Theta
f5 = figure();
for i = 1:length(theta_hat_toPlot(1,:))
    title_text = "θ_%d";
    subplot(ceil(length(theta_hat_toPlot(1,:)))/2,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'DisplayName','Predicted')
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , 'DisplayName','Real')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,main_folder + sub_name + "-theta" + '.jpeg')
close all