%% Offline design - Compute c_s, c_o, and alpha using P and K
% Initialize MATLAB interface
% clear;
close all;
% clc;

% Load P and K
load P_K_tracking.mat;
P(abs(P) < 1e-6) = 0; % round numerical tolerance small numbers to zero


%% Compute c_js for system constraints
% Obtain system (state and input) constraints
[~,~,~,~,L_s,~] = model.get_constraints();
L_s = L_s(2:2:end,:);
n_s = size(L_s,1);

% Compute constant factor used in c_s
c_s_const = inv(sqrtm(P))*[K',eye(model.n_x)];

% Compute c_s
c_s = zeros(n_s,1);
for i=1:n_s
  c_s(i) = norm(c_s_const*L_s(i,:)');
end


%% Compute c_o for obstacle avoidance constraints
% NOTE: can directly compute the 2-norm for matrix, since:
% - c_o only depends on L_o, not on l_o, and L_o can be scaled so ||L||=1
% - ||P^{-1/2} C' L'|| \leq ||P^{-1/2} C'|| ||L'|| \leq ||P^{-1/2} C'||
C = model.get_C();
c_o_const = inv(sqrtm(P))*C';
c_o = norm(c_o_const,2);


%% Determine terminal set scaling alpha
% Choose minimum safety distance d to obstacles
% Tune this to the smallest possible value d_min such that:
% (a) the TMPC is able to repeatedly find a solution ending in the terminal
% set
% (b) the PMPC can plan non-conservative collision-free trajectories
d = 0.1;  % 10 cm

% Compute alpha based on obstacle avoidance constraints
alpha = d/c_o;


%% Print results
% Determine how much space of half the constraint region is over after
% subtracting c_s*alpha
% Separate into lower and upper bounds of system constraints
sys_con_names = model.get_ux_names();
sys_con_lb = model.get_lower_bounds();
sys_con_ub = model.get_upper_bounds();

% Determine gaps between lower and upper bounds
sys_con_m = (sys_con_ub+sys_con_lb)/2;
sys_con_halfs = sys_con_ub-sys_con_m;

% Determine left overs of half the constraint regions
sys_con_halfs_left_overs = sys_con_halfs-c_s*alpha;

% Print results
fprintf("alpha: %3.5f\n\n", alpha);
fprintf("%s\n", "Left-over of half of the constraint region for:");
for i=1:size(sys_con_names,1)
  fprintf("%s:\t%3.3f / %3.3f\n",sys_con_names(i),sys_con_halfs_left_overs(i),sys_con_halfs(i));
end


%% Save results
save('offline_comps_tracking.mat','model','Q','R', ...
     'X','Y','eps_s','P','K','c_s', ...
     'timing', ...
     'c_o','d','alpha','sys_con_halfs_left_overs');

