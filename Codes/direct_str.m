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
%% direct with zero cancellation 
syms q;
main_folder_name = "1-1-direct-str";
main_folder = 'images/q1-final/' + main_folder_name + "/";
sub_name = "colored-noise-";
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
noise_poly = [1,0.4,0.2].';
deg_noise = length(noise_poly);
noise_variance = 0.005; % system noise variance
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

% desired B
Bplus = [BmPrime_main/BmPrime_main(1)];
Bminus = BmPrime_main(1);
Bm = conv(Bplus, Bminus);
BmPrime = BmPrime_main/BmPrime_main(1);


% initial R S T and Ao
R_estimated = [1,-0.1]; 
S_estimated = [1,0.1];
T_estimated = [1.1,0];
Ao = [1];

% R S and T which were calculated in the last part
R_real = [1,-0.3223,-0.4214];
S_real = [4.9308,-5.0634,1.5007];
T_real = [0.8293,0.2863,-0.002];


d0 = deg_desA - deg_desB;
[deg_R, deg_Ao] = find_degrees(Am, A, Ao, Bminus, Bplus, BmPrime);
deg_S = deg_R;
deg_T = deg_R; 

f_filter = conv(Am, Ao);

assert(f_filter(1)==1,"AmAo is not monic")

% initial conditions for y
skip_instances = max(deg_desA, deg_desB);

y = [];
y(1:skip_instances) = 0;
ym(1:skip_instances) = 0;
u(1:skip_instances) = 0;
uf(1:skip_instances) = 0;
yf(1:skip_instances) = 0;
ucf(1:skip_instances) = 0;

total_parameters = deg_R+1 + deg_S+1 + deg_T+1;


R_estimated_toPlot = zeros([num_samples, length(R_estimated)]);
S_estimated_toPlot = zeros([num_samples, length(S_estimated)]);
T_toPlot = zeros([num_samples, length(T_estimated)]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

theta_real = [A(2:end).'; B.'];
theta_desired = [Am(2:end).'; Bm.'];


theta_epsilon_zero = [0.5,0.5].';
els_solver = ELSClass(100 * eye(total_parameters+ length(theta_epsilon_zero)), 0.01 * ones([total_parameters,1]), theta_epsilon_zero);
for i = skip_instances:length(uc)
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].'; % assumed deg_desB = deg_B
    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly; %only for simulation. not involved in controller calculations directly
    
    y(i) = phi_t.' * theta_real + noise_t + B * v(i:-1:i-(length(B)-1)).';
    u(i) = T_estimated * [uc(i:-1:i-(length(T_estimated)-1))].' + S_estimated * [-y(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [u(i-1:-1:i-(length(R_estimated)-1))].';
    u(i) = u(i)./R_estimated(1);

    uf(i) = u(i) - f_filter(2:end) * uf(i-1:-1:i-(length(f_filter)-1)).';
    yf(i) = y(i) - f_filter(2:end) * yf(i-1:-1:i-(length(f_filter)-1)).';
    ucf(i) = uc(i) - f_filter(2:end) * ucf(i-1:-1:i-(length(f_filter)-1)).';

    phi_d0_filtered = [uf(i-d0:-1:i-d0-deg_R), yf(i-d0:-1:i-d0-deg_S), -ucf(i-d0:-1:i-d0-deg_T)].';
    
    phi_t_m = [-y(i-1:-1:i-(deg_desA - 1)), uc(i-1:-1:i-deg_desB)].';
    ym(i) = phi_t_m.' * theta_desired;
    error = y(i) - ym(i);

    theta_hat_new = els_solver.update_ELS(error, phi_d0_filtered);
   
    R_estimated = theta_hat_new(1:deg_R+1).';
    S_estimated = theta_hat_new(deg_R+2:deg_R + deg_S+2).';
    T_estimated = theta_hat_new(deg_R + deg_S+3:deg_R + deg_S+deg_T+3).';


    R_estimated_toPlot(i, :) = R_estimated;
    S_estimated_toPlot(i, :) = S_estimated;
    T_toPlot(i, :) = T_estimated;
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

% R 
% total_plot_rows = ceil(length(R_solved)/2) + ceil(length(S_solved)/2)+ ceil(length(T)/2) - 1;
% f2 = figure();
% f2.Position = [0 0 400 900];
% for i = 2:(length(R_solved)) % first parameter of R is 1
%     subplot(total_plot_rows,2,i-1);
%     title_text = "R_%d";
%     plot(1:num_samples, R_solved_toPlot(:,i), 'DisplayName','Predicted') 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
%     title(sprintf(title_text, i-1));
%     legend('Location','best');
%     xlabel("sample number")
% 
% end
% end_of_r_plot = i-1;
% 
% % S
% for i = 1:(length(S_solved)) % first parameter of R is 1
%     subplot(total_plot_rows,2,i + end_of_r_plot);
%     title_text = "S_%d";
%     plot(1:num_samples, S_solved_toPlot(:,i), 'DisplayName','Predicted') 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
%     title(sprintf(title_text, i));
%     legend('Location','best');
%     xlabel("sample number")
% end
% 
% end_of_S_plot = i + end_of_r_plot;
% % T
% for i = 1:(length(T)) % first parameter of R is 1
%     subplot(total_plot_rows,2,i+end_of_S_plot);
%     title_text = "T_%d";
%     plot(1:num_samples, T_toPlot(:,i), 'DisplayName','Predicted') 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * T_real(i), 'DisplayName','Real') 
%     title(sprintf(title_text, i));
%     legend('Location','best');
%     xlabel("sample number")
% end
% saveas(gcf,'images/q2/' + main_title+ "-T" + '.jpeg')
% close all

% Theta
% f5 = figure();
% for i = 1:length(theta_hat_toPlot(1,:))
%     title_text = "θ_%d";
%     subplot(ceil(length(theta_hat_toPlot(1,:)))/2,2,i);
%     plot(1:num_samples, theta_hat_toPlot(:,i), 'DisplayName','Predicted')
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , 'DisplayName','Real')
%     title(sprintf(title_text, i));
%     legend('Location','best');
%     xlabel("sample number")
% end
% saveas(gcf,main_folder + sub_name + "theta" + '.jpeg')
% close all
disp("finished")

