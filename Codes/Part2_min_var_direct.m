clear;

Ts = 0.2;
sysD = tf([1,-0.8],[1,-1.5,0.6],Ts);
[numD, denD] = tfdata(sysD, 'v');


% bodeplot(sysD);

B = numD;
A = denD;
% step(sysD)
%%
if B(1) == 0
    B = B(2:end);
end



%% direct with zero cancellation 
rng(50);

main_folder_name = "2-1-direct-minVar";
main_folder = 'images/q1-final/' + main_folder_name + "/";
sub_name = "direct";
main_title = sub_name;
close_all = true;

% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;

% noise and its degree
% noise_poly = [1,0.4,0.2].';
noise_poly = [1,-0.8,0.1].';
len_noise = length(noise_poly);
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
d0 = deg_desA - deg_desB;
A_estimated = [1,0,0]; % initial conditions
B_estimated = [0.01, 0.035];
% RLS initial parameters (C)
theta_epsilon_zero = [-0.5,0.5].';
C_estimated = [1, theta_epsilon_zero.'];

assert(length(A_estimated)==len_desA && length(B_estimated)==len_desB, "initial condistions are not true")


% initial R S T and Ao
R_estimated = [1,-0.1]; 
S_estimated = [1, -0.1];
len_R = length(R_estimated);
len_S = length(S_estimated);

% R S and T which were calculated in the last part
R_real = [1.0000   -0.8000];
S_real = [0.7000   -0.5000];
% T_real = [0.8293,0.2863,-0.002];
P_filter = noise_poly.';
Q_filter = [1];


% initial conditions for y
skip_instances = max(len_desA, len_desB);

y = [];
y(1:skip_instances) = 0;
ym(1:skip_instances) = 0;
u(1:skip_instances) = 0;
uf(1:skip_instances) = 0;
yf(1:skip_instances) = 0;


total_parameters = len_R + len_S;


R_estimated_toPlot = zeros([num_samples, length(R_estimated)]);
S_estimated_toPlot = zeros([num_samples, length(S_estimated)]);

theta_hat_toPlot = zeros([num_samples, total_parameters]);

theta_real = [A(2:end).'; B.'];


els_solver = ELSClass(100 * eye(total_parameters+ length(theta_epsilon_zero)), 0.01 * ones([total_parameters,1]), theta_epsilon_zero);
for i = skip_instances:num_samples
    phi_t = [-y(i-1:-1:i-(len_desA - 1)), u(i-1:-1:i-len_desB)].';  

    noise_t = [noise(i:-1:i-(len_noise-1))] * noise_poly;
    y(i) = phi_t.' * theta_real + noise_t + B * v(i:-1:i-(length(B)-1)).';

    u(i) = (S_estimated * [-y(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [u(i-1:-1:i-(length(R_estimated)-1))].')/R_estimated(1);

    uf(i) = Q_filter * u(i:-1:i-length(Q_filter)+1).' - P_filter(2:end) * uf(i-1:-1:i-(length(P_filter)-1)).';
    yf(i) = Q_filter * y(i:-1:i-length(Q_filter)+1).' - P_filter(2:end) * yf(i-1:-1:i-(length(P_filter)-1)).';

    phi_filtered = [uf(i-d0:-1:i-d0-len_R+1), yf(i-d0:-1:i-d0-len_S+1)].';
    ym(i) = 0;
    error = y(i) - ym(i);

    theta_hat_new = els_solver.update_ELS(error, phi_filtered);
   
    R_estimated = theta_hat_new(1:len_R).';
    S_estimated = theta_hat_new(len_R+1:len_R + len_S).';


    R_estimated_toPlot(i, :) = R_estimated;
    S_estimated_toPlot(i, :) = S_estimated;
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
total_plot_rows = ceil(length(R_estimated)/2) + ceil(length(S_estimated)/2);
f2 = figure();
for i = 1:(length(R_estimated)) % first parameter of R is 1
    subplot(total_plot_rows,2,i);
    title_text = "R_%d";
    plot(1:num_samples, R_estimated_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")

end
end_of_r_plot = i;

% S
for i = 1:(length(S_estimated)) % first parameter of R is 1
    subplot(total_plot_rows,2,i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_estimated_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,main_folder + sub_name + "-S-R" + '.jpeg')
close all
