classdef HovergamesModel
  % SYSTEM Hovergames model class
  %
  %   Args:
  %     None
  %
  %   Properties:
  %     n_u : Number of inputs
  %     n_x : Number of states
  %     t_s : Sampling time
  %
  %   Inputs:
  %     phi_c    : Commanded roll angle (using ZYX Euler angles)
  %     theta_c  : Commanded pitch angle (using ZYX Euler angles)
  %     psi_c    : Commanded yaw angle (using ZYX Euler angles)
  %     thrust_c : Commanded thrust (thrust_c = g for hovering)
  %
  %   States:
  %     x       : Position along x-axis (world frame)
  %     y       : Position along y-axis (world frame)
  %     z       : Position along z-axis (world frame)
  %     vx      : Velocity along x-axis (world frame)
  %     vy      : Velocity along y-axis (world frame)
  %     vz      : Velocity along z-axis (world frame)
  %     phi     : Roll angle (using ZYX Euler angles)
  %     theta   : Pitch angle (using ZYX Euler angles)
  %     psi     : Yaw angle (using ZYX Euler angles)
  %     thrust  : Thrust (body frame, thrust = g for hovering)
  %
  %   System constraints:
  %     L_u : Input constraint weights (L_u * u <= l_u)
  %     l_u : Input constraint bounds (L_u * u <= l_u)
  %     L_x : State constraint weights (L_x * x <= l_x)
  %     l_x : State constraint bounds (L_x * x <= l_x)
  %     L   : Constraint weights (L = [L_u;L_x]
  %     l   : Constraint bounds (l = [l_u;l_x])
  %

  properties (Constant)
    % Note: all quantities are given in SI units
    % Model dimensions
    n_x = 10;
    n_u = 4;

    % Constants
    g = 9.81;

    % Attitude and thrust models fitted on real drone
    A_phi        = -5.55;
    B_phi        = 5.55;
    A_theta      = HovergamesModel.A_phi;
    B_theta      = HovergamesModel.B_phi;
    A_psi        = -1.773;
    B_psi        = 1.773;
    A_thrust     = -20;
    B_thrust     = 20;

    % Position selection matrix
    n_pos = 2;
    C = [eye(HovergamesModel.n_pos),zeros(HovergamesModel.n_pos,HovergamesModel.n_x-HovergamesModel.n_pos)];

    % System constraints
    att_min = -pi/6;
    att_max = pi/6;
    thrust_min = 5;
    thrust_max = 15;

    % Input constraints
    L_u = kron(eye(HovergamesModel.n_u),[-1;1]);
    l_u = [-HovergamesModel.att_min;HovergamesModel.att_max;
           -HovergamesModel.att_min;HovergamesModel.att_max;
           -HovergamesModel.att_min;HovergamesModel.att_max;
           -HovergamesModel.thrust_min;HovergamesModel.thrust_max];

    % State constraints
    L_x = kron(eye(HovergamesModel.n_x),[-1;1]);
    l_x = [15;15;
           15;15;
           0;4;
           2;2;
           2;2;
           2;2;
           -HovergamesModel.att_min;HovergamesModel.att_max;
           -HovergamesModel.att_min;HovergamesModel.att_max;
           -HovergamesModel.att_min;HovergamesModel.att_max;
           -HovergamesModel.thrust_min;HovergamesModel.thrust_max];

    % Grid settings
    % n_grid : default number of grid points
    % n_grid_vert : number of grid points for taking vertices
    % To construct LMIs
    n_grid_construct = 5;
    n_grid_vertices_construct = 2;
    % To check whether LMIs hold at a finer grid in the input-state space
    n_grid_check = 21;
    n_grid_vertices_check = 2;
  end

  methods(Static)
    function names = get_ux_names()
      names = ["phi_c";
               "theta_c";
               "psi_c";
               "thrust_c";
               "x";
               "y";
               "z";
               "vx";
               "vy";
               "vz";
               "phi";
               "theta";
               "psi";
               "thrust"];
    end
  end

  methods
    % Continous-time model
    function f = fun(obj,x,u)
      % Assign current states
      vx     = x(4);
      vy     = x(5);
      vz     = x(6);
      phi    = x(7);
      theta  = x(8);
      psi    = x(9);
      thrust = x(10);

      % Assign inputs
      phi_c    = u(1);
      theta_c  = u(2);
      psi_c    = u(3);
      thrust_c = u(4);

      % Calculate system update
      f = [vx;...
           vy;...
           vz;...
           thrust*(sin(psi)*sin(phi)+cos(phi)*sin(theta)*cos(psi));...
           thrust*(-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi));...
           thrust*(cos(phi)*cos(theta))-obj.g;...
           obj.A_phi*phi+obj.B_phi*phi_c;...
           obj.A_theta*theta+obj.B_theta*theta_c;...
           obj.A_psi*psi+obj.B_psi*psi_c;...
           obj.A_thrust*thrust+obj.B_thrust*thrust_c;...
          ];
    end

    function C = get_C(obj)
      C = obj.C;
    end

    % Get system constraints in matrix form
    function [L_u,l_u,L_x,l_x,L,l] = get_constraints(obj)
      L_u = obj.L_u;
      l_u = obj.l_u;
      L_x = obj.L_x;
      l_x = obj.l_x;

      % Combined state and input constraints
      L = [L_u,zeros(2*obj.n_u,obj.n_x);
           zeros(2*obj.n_x,obj.n_u),L_x];
      l = [l_u;l_x];
    end

    % Get casadi derivative functions
    % Note: use full(A_fun(x,u)) to obtain a normally-sized A matrix
    function [f_fun,A_fun,B_fun] = get_fABfun(obj)
      import casadi.*
      x_sym = MX.sym('x',obj.n_x,1);
      u_sym = MX.sym('u',obj.n_u,1);
      f     = obj.fun(x_sym,u_sym);
      f_fun = casadi.Function('f_fun',{x_sym,u_sym},{f});
      A     = jacobian(f,x_sym); 
      B     = jacobian(f,u_sym); 
      A_fun = casadi.Function('A_fun',{x_sym,u_sym},{A});
      B_fun = casadi.Function('B_fun',{x_sym,u_sym},{B});
    end

    function lb = get_lower_bounds(obj)
      lb = [-obj.l_u(1:2:end);
            -obj.l_x(1:2:end)];
    end

    function ub = get_upper_bounds(obj)
      ub = [obj.l_u(2:2:end);
            obj.l_x(2:2:end)];
    end

    % Compute Jacobian on grid
    % This is tailored to specific nonlinear dynamics and needs to be
    % adjusted to the model implemented in this class
    function [ugrid,xgrid,Agrid,Bgrid,n] = get_uxABgrid(obj,do_print,do_check)
      n_bytes = 0;
      n_grid = obj.n_grid_construct;
      n_grid_vertices = obj.n_grid_vertices_construct;
      if do_check
        n_grid = obj.n_grid_check;
        n_grid_vertices = obj.n_grid_vertices_check;
      end
      % Gridding choices:
      % - Thrust enters affinely in f -> linearly in jacobian -> take vertices
      % - Angles: psi,phi,theta: very nonlinear -> grid
      % - Other states and inputs do not enter the Jacobians
      n = n_grid^3*n_grid_vertices;
      ugrid = zeros(obj.n_u,n);
      xgrid = zeros(obj.n_x,n);
      Agrid = zeros(obj.n_x,obj.n_x,n);
      Bgrid = zeros(obj.n_x,obj.n_u,n);
      [~,A_fun,B_fun] = obj.get_fABfun();
      i = 0;
      for phi = linspace(HovergamesModel.att_min,HovergamesModel.att_max,n_grid)
        for theta = linspace(HovergamesModel.att_min,HovergamesModel.att_max,n_grid)
          for psi = linspace(HovergamesModel.att_min,HovergamesModel.att_max,n_grid)
            for thrust = linspace(HovergamesModel.thrust_min,HovergamesModel.thrust_max,n_grid_vertices)
              i = i + 1;
              u = zeros(obj.n_u,1);
              x = [zeros(6,1);phi;theta;psi;thrust];
              ugrid(:,i) = u;
              xgrid(:,i) = x;
              Agrid(:,:,i) = full(A_fun(x,u));
              Bgrid(:,:,i) = full(B_fun(x,u));
              if do_print && rem(i,10) == 0
                fprintf(repmat('\b',1,n_bytes));
                n_bytes = fprintf("Constructed gridding point %i/%i\n",i,n);
              end
            end
          end
        end
      end
      if do_print
        fprintf(repmat('\b',1,n_bytes));
        fprintf("Constructed gridding point %i/%i\n",i,n);
      end
    end
  end
end
