% Define parameters
T = 30 * 60;
w = 2 * pi / T; % Angular velocity of London Eye in rad/s
r = 60;  % Radius of London Eye in meters
g = 9.81;     % Gravitational acceleration in m/s^2

% Define a range of theta (angle in radians)

% Define the angle theta as a function of time
t = linspace(0, 30, 100);
theta = w * t * 60 - pi/2;

% Compute g-force
Gx = - 4 * pi * r * cos(theta) / (g * T^2);
Gz = 1 - 4 * pi * r * sin(theta) / (g * T^2);

% Calculate the maximum and minimum g-force
Gx_max = max(Gx);
Gx_min = min(Gx);
Gz_max = max(Gz);
Gz_min = min(Gz);

% Display the maximum and minimum g-force in the command window
disp(['Maximum Gx: ', num2str(Gx_max)]);
disp(['Minimum Gx: ', num2str(Gx_min)]);
disp(['Maximum Gz: ', num2str(Gz_max)]);
disp(['Minimum Gz: ', num2str(Gz_min)]);

% Subplot for Gx
subplot(2, 1, 1);
plot(t, Gx, 'LineWidth', 2); 
set(gca, 'FontSize', 15);
xlabel('t (min)', 'FontSize', 18);
ylabel('$G_x$', 'Interpreter', 'latex', 'FontSize', 18);
legend('$G_x$', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
axis tight;

% Subplot for Gz
subplot(2, 1, 2);
plot(t, Gz, 'r', 'LineWidth', 2); 
set(gca, 'FontSize', 15);
xlabel('t (min)', 'FontSize', 18);
ylabel('$G_z$', 'Interpreter', 'latex', 'FontSize', 18);
legend('$G_z$', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
axis tight;