%% chebivp.m -- an executable m-file for solving an initial-value problem
% Automatically created in CHEBGUI by user maye.

%% Problem description.
% Solving
%   u'' + sin(u) - cos(u)*sin(u+4/pi) = 0,
% for t in [0 10], subject to
%   u'(0) = 0,
%   u(0) = 0.

%% Problem set-up.
% Define the domain.
dom = [0 10];

M = 500;
l = 1.5;
v_max = 20;
angle = 0 * pi / 180;
A0 = 4*l.^2;
g = 9.81;
Cd = 1.12;
ro = 1.29;
t1 = 0.5;
t2 = 1.5;
t3 = 3;

vw = chebfun({...
        @(t) (v_max / t1) * t, ...             % From 0 to t1: Linear increase
        @(t) v_max, ...                         % From t1 to t2: Stay v_max
        @(t) v_max * (1 - (t - t2) / (t3 - t2)), ... % From t2 to t3: Linear decrease
        @(t) 0}, ...                           % From t3 to end: Zero velocity
        [0, t1, t2, t3, dom(2)]);              

% Assign the differential equation to a chebop on that domain.
N = chebop(@(t, u) M * l * diff(u, 2) + M * g * sin(u) - 0.5 * sqrt(2) * A0 * Cd * ro * (vw.^2) .* cos(u-angle) .* sin((u + 4 ./ pi)), dom);

% Set up the rhs of the differential equation so that N(u) = rhs.
rhs = 0;

% Assign boundary conditions to the chebop.
N.lbc = @(u) [diff(u); u]; % Corrected to explicitly set boundary conditions.

%% Setup preferences for solving the problem.
% Create a CHEBOPPREF object for passing preferences.
% (See 'help cheboppref' for more possible options.)
options = cheboppref();

% Specify the IVP solver to use. Possible options are:
%   Time-stepping solvers:
%     'ode113' (default), 'ode15s' or 'ode45'.
%   Global methods:
%     'values' or 'coefficients'.
options.ivpSolver = 'ode113';

%% Solve!
% Call solveivp() to solve the problem.
% (With the default options, this is equivalent to u = N\rhs.)
u = solveivp(N, rhs, options); % u repersent angle beta
w = diff(u); % angular velocity
diff_w = diff(w); 

% Calculate gforce using element-wise operations
ar = l * (w).^2;
at = l * diff_w;
Gx = sin(u) + at/g;
Gz = cos(u) + ar/g;

% Calculate the maximum and minimum g-force
Gx_max = max(Gx);
Gx_min = min(Gx);
Gz_max = max(Gz);
Gz_min = min(Gz);

% Calculate the maximum velocity
v = w .* l;
vmax = max(v);

% Display the maximum and minimum g-force in the command window
disp(['Maximum Gx: ', num2str(Gx_max)]);
disp(['Minimum Gx: ', num2str(Gx_min)]);
disp(['Maximum Gz: ', num2str(Gz_max)]);
disp(['Minimum Gz: ', num2str(Gz_min)]);
disp(['Velocity: ', num2str(vmax)]);

% Check if the maximum or minimum value of u exceeds the specified limits
u_max = max(u);
u_min = min(u);
if u_max > pi/3 || u_min < -pi/3
    disp(['beta: ', num2str(u_max)]);
    fprintf('Stopping due to u exceeding limits for beta');
end

%% Plotting the results

% Figure for angle and angle velocity over time
figure;
plot(u, 'LineWidth', 2); 
hold on;
plot(w, 'LineWidth', 2); 
hold off;
% Set axis properties
set(gca, 'FontSize', 15);
grid on;
% Labels
xlabel('t', 'FontSize', 18);
ylabel('$\beta\ \&\ w$', 'Interpreter', 'latex', 'FontSize', 18);
% Legend
legend({'$\beta$', '$w$'}, 'Interpreter', 'latex', 'FontSize', 18);

% Figure for g-force over time
figure;
plot(Gz, 'LineWidth', 2);
% Set axis properties
set(gca, 'FontSize', 15);
grid on;
% Labels
xlabel('t', 'FontSize', 18);
ylabel('$G_z$', 'Interpreter', 'latex', 'FontSize', 18);
% Legend
legend('$G_z$', 'Interpreter', 'latex', 'FontSize', 18);

figure;
plot(Gx, 'LineWidth', 2);
% Set axis properties
set(gca, 'FontSize', 15);
grid on;
% Labels
xlabel('t', 'FontSize', 18);
ylabel('$G_x$', 'Interpreter', 'latex', 'FontSize', 18);
% legend
legend('$G_x$', 'Interpreter', 'latex', 'FontSize', 18);