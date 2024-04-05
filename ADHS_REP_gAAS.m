% ADHS_REP_gAAS: MATLAB script for computing a gAAS relation for the case study of the vehicle, simulating the system, and plotting figure 3 and 4 in the paper
% Code for the Paper entitled "Hierarchical Control for Cyber-Physical Systems via General Approximate Alternating Simulation Relations" in ADHS 2024
%   Authors: Bingzhuo Zhong
%   Date: April 1, 2024

clear all

% add necessary tools for the following computation
addpath('Tools');

%% Part 1 Computing a gAAS relation for the case study of the vehicle 

% instantiating an object for computation of gAAS relation
car = gAAS_relation_options();

% specifying system dynamics for the concrete and abstract systems
dt = 0.1;
car.A = [1 dt;0 1];
car.B = [dt^2/2;dt];
car.C = [1 0];
car.A_hat = 1;
car.B_hat = dt;
car.C_hat = 1;

% defining the state and input sets of the abstract system
X_hat_temp.A = [1;-1];
X_hat_temp.b = [1000;1000];
car.X_hat = {X_hat_temp};   % Cell containing serveral polytopes

U_hat_temp.A = [1;-1];
U_hat_temp.b = [35;35];
car.U_hat = {U_hat_temp};   % Cell containing serveral polytopes

% specifying matrices P, Q, and D in equation (8)
car.P = [1;0];
car.Q = 0;
car.D = [0;1];

% varepsilon in condition (ii) in Definition 3
car.varepsilon = 0.02;

% solving the gAAS relation according to equations (13) and (14)
kappa = 0.65;

% for repeatability, fix K and compute M in Lemma 4.1, then one can get the same results reported in the paper
K_fix = [-264 -19.3];
[fea,car] = car.solve_kappa_fixK(kappa,K_fix);

% one can also try to use the following code to compute M and K in Lemma 4.1 jointly
% [fea,car] = car.solve_kappa(kappa);

% save the gAAS relation in a .mat data file
save('ADHS_2D_Result','car')

%% Part 2: Simulating the case study

% computing the dimensions of the state and input sets for the concrete and abstract systems
x_dim = size(car.A,2);
xhat_dim = size(car.A_hat,2);
u_dim = size(car.B,2);
uhat_dim = size(car.B_hat,2);
y_dim = size(car.C,2);

% setting the simulation time horizon
time_horizon = 500;

% picking initial state of the concrete and the abstract system
init_x0 = [0;0];
init_xhat0 = 0;
init_uhat = 0;

% specifying K for the high level controller
K_h = 10;

% initializing matrces for storaging the simulation results
x_sequence = zeros(x_dim,time_horizon+1);               % sequance of states of the concrete system
xhat_sequence = zeros(xhat_dim,time_horizon+1);         % sequance of states of the abstract system
u_sequence = zeros(u_dim,time_horizon+1);               % sequance of control inputs for the concrete system
uhat_sequence = zeros(uhat_dim,time_horizon);           % sequance of control inputs for the abstract system
y_sequence = zeros(y_dim,time_horizon+1);               % sequance of outputs of the concrete system
yhat_sequence = zeros(y_dim,time_horizon+1);            % sequance of outputs of the abstract system
e_sequence = zeros(1,time_horizon+1);                   % sequance of errors between the output of the concrete and the abstract systems

% initializing the system
x_sequence(:,1) = init_x0;
xhat_sequence(:,1) = init_xhat0;
u_hat_cur = init_uhat;

% get a series of set-points for the abstract system by calling the function trajectory_gen.m in the file Tools
x_series = trajectory_gen(car,time_horizon+1,3);
x_hat_set = x_series(1,:);

for i=1:1:time_horizon
    % simulation
    
    % initialization
    x_cur = x_sequence(:,i);
    x_hat_cur = xhat_sequence(:,i);
    y_sequence(:,i) = car.C*x_cur;
    yhat_sequence(:,i) = car.C_hat*x_hat_cur;
    e_sequence(:,i) = sqrt((y_sequence(:,i)-yhat_sequence(:,i))'*(y_sequence(:,i)-yhat_sequence(:,i)));

    % get the controller for the abstract system
    if i ==1
        u_hat_cur = init_uhat;
    else
        u_sug = -K_h*(x_hat_cur-x_hat_set(:,i));

        % get the controller considering the changing rate of the controller for the abstract system as discussed in Lemma 4.1
        u_hat_cur = car.provide_u(uhat_sequence(:,i-1),u_sug);
    end

    % recording the control input for the abstract system
    uhat_sequence(:,i) = u_hat_cur;

    % computing the refinement controller for the concrete system
    u_cur = car.K*(x_cur-car.P*x_hat_cur-car.S*u_hat_cur)+car.Q*x_hat_cur+car.R*u_hat_cur;
    u_sequence(:,i) = u_cur;

    % update the state of the concrete system
    x_nxt= car.A*x_cur+car.B*u_cur;

    % update the state of the abstract system
    x_hat_nxt= car.A_hat*x_hat_cur+car.B_hat*u_hat_cur;

    % save the next state
    x_sequence(:,i+1) = x_nxt;
    xhat_sequence(:,i+1) = x_hat_nxt;

end
% record the output and error at the last step
y_sequence(:,time_horizon+1) = car.C*x_nxt;
yhat_sequence(:,time_horizon+1) = car.C_hat*x_hat_nxt;
e_sequence(:,time_horizon+1) = sqrt((y_sequence(:,time_horizon+1)-yhat_sequence(:,time_horizon+1))'*(y_sequence(:,time_horizon+1)-yhat_sequence(:,time_horizon+1)));

% save the simulation results in a .mat data file
save('ADHS_sim_results','time_horizon','y_sequence','yhat_sequence','e_sequence','uhat_sequence','u_sequence','x_sequence','x_series','xhat_sequence')

%% Part 3: Plotting Figure 3 and 4 in the paper according to the simulation result

figure(4)
% plotting the evolution of y and yhat in Figure 4
subplot(2,1,1)
t1 = 0:1:time_horizon+1;
data1 = y_sequence(1,:);
[tn,xn] = pwc1_plot(t1,data1);
plot(tn,xn,'b','LineWidth',2);
hold on
t2 = 0:1:time_horizon+1;
data2 = yhat_sequence(1,:);
[tn2,xn2] = pwc1_plot(t2,data2);
plot(tn2,xn2,'r--','LineWidth',2);
legend('$y(k)$','$\hat{y}(k)$','Location','northeast','Interpreter','Latex','Fontsize',30);
xlabel('$k$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
ylabel('$y\backslash \hat{y}$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
set(gca,'linewidth',2,'fontsize',40,'fontname','Times');
axis([0 500 0 Inf])

% plotting the evolution of the error between y and yhat in Figure 4
subplot(2,1,2)
t3 = 0:1:time_horizon+1;
data3 = e_sequence(1,:);
[tn3,xn3] = pwc1_plot(t3,data3);
plot(tn3,xn3,'b','LineWidth',2);
xlabel('$k$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
ylabel('$|| y- \hat{y}||$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
set(gca,'linewidth',2,'fontsize',40,'fontname','Times');
axis([0 500 -0.0005 0.001])


figure(3)
% plotting evolution of the nuhat in Figure 3
subplot(2,1,1)
t4 = 0:1:time_horizon;
data4 = uhat_sequence(1,:);
[tn4,xn4] = pwc1_plot(t4,data4);
plot(tn4,xn4,'b','LineWidth',2);
xlabel('$k$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
ylabel('$\hat{\nu}$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
set(gca,'linewidth',2,'fontsize',40,'fontname','Times');
axis([0 499 -35 35])
xticks([0:50:450 499])
yticks([-35 -15 0 15 35])

% plotting evolution of the change of the nuhat in Figure 3
subplot(2,1,2)
t5 = 0:1:time_horizon-1;
data5 = uhat_sequence(1,2:time_horizon)-uhat_sequence(1,1:time_horizon-1);
[tn5,xn5] = pwc1_plot(t5,data5);
plot(tn5,xn5,'b','LineWidth',2);
xlabel('$k$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
ylabel('$\Delta \hat{\nu}$','Interpreter','Latex','Fontsize',60,'FontWeight','bold','fontname','Italic');
set(gca,'linewidth',2,'fontsize',40,'fontname','Times');
axis([0 498 -0.36 0.36])
xticks([0:50:450 498])
yticks(-0.36:0.18:0.36)