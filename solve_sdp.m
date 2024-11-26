%% Offline design - solve SDP to obtain X, Y, and c_s
% Compute quadratic incremental Lyapunov function using SDP

% Initialize MATLAB interface
clear;
close all;
clc;

% Import necessary libraries
import casadi.* % CasADi is used to automatically compute derivates of nonlinear dynamic

% Start timing
timer = tic;

% Select model
model = HovergamesModel();

% Check solution (boolean actives additional validation of solution after
% solving SDP with tolerance solver_tol and finer grid)
do_check_sol = true;
solver_tol = 1e-6;

% Print information
do_print = true;


%% Set controller information
% NOTE:
% - For simulation with model mismatch, tune these values the following
% way:
% - Increase Q position weights to track positions despite model mismatch
% - Increase R attitude command weights to prevent the system from becoming
% unstable if the model is not perfect. The thrust model is less accurate
% and more critical, so do not increase this value too much
Q = diag([100,100,100,1,1,1,5,5,5,5]);
R = diag([100,100,100,5]);

Q_eps = sqrtm(Q);
R_eps = sqrtm(R);


%% Describe system constraints sets
% System (state and input) constraints as polytopic constraint
[~,~,~,~,L_s,~] = model.get_constraints();
L_s = L_s(2:2:end,:); % this uses the fact that our constraints are always double sided
n_s = size(L_s,1);

% Separate into lower and upper bounds of system constraints
sys_con_lb = model.get_lower_bounds();
sys_con_ub = model.get_upper_bounds();

% Determine gaps between lower and upper bounds
% NOTE: this is for scaling eps_s=c_s^2 in the cost, since otherwise the
% terminal set becomes very small for some inputs/states with respect to 
% other inputs/states
sys_con_m = (sys_con_ub+sys_con_lb)/2;
sys_con_halfs = sys_con_ub-sys_con_m;


%% Set up optimization problem
% Define decision variables
eps_s = sdpvar(n_s,1); % this variable corresponds to c_s^2
X = sdpvar(model.n_x);
Y = sdpvar(model.n_u,model.n_x);

% Define objective
% NOTE:
% - Infeasible for (c <= 0 (unbounded), c >= 1e7 (numerical issues))
% - Penalize constraints tightening to descrease terminal set size
% - Normalize constraints w.r.t. their constraints interval to tighten each
% constraint equally (as a percentage of the complete range)
% - c large => reshape P to achieve equal tightening in every dimension
% - Tune c, starting from a low value until the input tightening saturates
c = 20;
obj = -log(det(X)) + c*sum(eps_s./sys_con_halfs.^2);

% Define solver options
ops = sdpsettings('solver','sdpt3','verbose',1);

timing.t_init = toc(timer);
if do_print
  fprintf("Init time: %f\n\n",timing.t_init);
end

% Define common constraint satisfaction LMIs
timer_lmi = tic;
if do_print
  fprintf("Constructing LMIs\n");
end
do_check = false; % this boolean changes which grid points are used (SDP vs. validation)
[~,~,Agrid,Bgrid,n_uxABgrid] = model.get_uxABgrid(do_print,do_check); % compute linearization
lmi_track = construct_check_tracking_LMIs(do_print,do_check,model,Q_eps,R_eps,Agrid,Bgrid,n_uxABgrid,X,Y,solver_tol); % LMI for terminal cost
lmi_sys = construct_check_sys_constraint_LMIs(do_print,do_check,n_s,L_s,eps_s,X,Y,solver_tol); % LMI for constraint tightening c_s
con = [lmi_track;lmi_sys];
timing.t_lmi = toc(timer_lmi);
if do_print
  fprintf("LMI computation time: %f\n\n",timing.t_lmi);
end


%% Optimize
timer_opt = tic;
if do_print
  fprintf("Optimize");
end
diagnostics = optimize(con,obj,ops);
timing.t_opt = toc(timer_opt);
if do_print
  print_diagnostics(diagnostics);
  fprintf("Optimization time: %f\n\n",timing.t_opt);
end


%% Obtain solution
timer_sol = tic;
eps_s = value(eps_s);
c_s = sqrt(eps_s);
X = value(X);
Y = value(Y);
P = inv(X);
K = Y/X;

if do_print
  fprintf("cond(P): %e\n", cond(P));
  fprintf("max eig(P): %e\n", max(eig(P),[],"all"));
end

timing.t_sol = toc(timer_sol);
if do_print
  fprintf("Solution time: %f\n\n",timing.t_sol);
end


%% Check solution
% Checks solution on finer grid and gives message if tolerance is exceeded
timer_check_sol = tic;
if do_check_sol
  do_check = true;
  if do_print
    fprintf("Checking LMIs\n");
  end
  [ugrid,xgrid,Agrid,Bgrid,n_uxABgrid] = model.get_uxABgrid(do_print,do_check_sol);
  construct_check_tracking_LMIs(do_print,do_check,model,Q_eps,R_eps,Agrid,Bgrid,n_uxABgrid,X,Y,solver_tol);
  construct_check_sys_constraint_LMIs(do_print,do_check,n_s,L_s,eps_s,X,Y,solver_tol);
end
timing.t_check_sol = toc(timer_check_sol);
if do_print
  fprintf("Solution checking time: %f\n\n",timing.t_check_sol);
end


%% Stop total timing
timing.t_total = toc(timer);
if do_print
  fprintf("Total time: %f\n\n",timing.t_total);
end


%% Store results
save('P_K_tracking.mat','model','Q','R', ...
     'X','Y','eps_s','P','K','c_s', ...
     'timing');

