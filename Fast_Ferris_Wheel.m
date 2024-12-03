%% Equation solve for Fast Ferris Wheels

% Define parameters
T = 6;                        % Period time
w_max = 2 * pi / T;       % Maximum angular velocity of Ferris Wheel in rad/s
r = 16;                          % Radius of Fast Ferris Wheel in meters
g = 9.81;                       % Gravitational acceleration in m/s^2

% Define the parameters of the piecewise function
t1 = 25;                         
t2 = 50;                        
t3 = 75;                        

% Define the piecewise function for angular velocity over time
f = @(t) (t >= 0 & t < t1) .* (w_max / t1 .* t) + ...                
         (t >= t1 & t < t2) .* w_max + ...                           
         (t >= t2 & t <= t3) .* (w_max - (w_max / (t3 - t2)) .* (t - t2));

% Define the piecewise function for derivative angular velocity over time
h = @(t) (t >= 0 & t < t1) .* (w_max / t1) + ...                
         (t >= t1 & t < t2) .* 0 + ...                           
         (t >= t2 & t <= t3) .* (- (w_max / (t3 - t2))); 

% Define the time range for calculating and plotting
t = linspace(0, t3, 1000);      

% Calculate the angular velocity and its derivatives at each time point
w = f(t);
diff_w = h(t);

% Define the angle theta as a function of time
F = @(t) (t >= 0 & t < t1) .* (w_max / (2 * t1) .* t.^2 - pi/2) + ...
         (t >= t1 & t < t2) .* (w_max .* t - (w_max * t1 / 2) - pi/2) + ...
         (t >= t2 & t <= t3) .* (w_max .* t - (w_max / (2 * (t3 - t2))) .* (t.^2 - 2 * t2 * t) ...
         - (w_max * t2^2) / (2 * (t3 - t2)) - (w_max * t1 / 2) - pi/2);

theta = F(t);

% Compute g-force over time
Gx = - (diff_w.^2 .* r .* sin(theta) + w.^2 .* r .* cos(theta)) / g;
Gz = 1 + (diff_w.^2 .* r .* cos(theta) - w.^2 .* r .* sin(theta)) / g;

% Calculate the maximum and minimum g-force
Gx_max = max(Gx);
Gx_min = min(Gx);
Gz_max = max(Gz);
Gz_min = min(Gz);

% Display the maximum anum2strnd minimum g-force in the command window
disp(['Maximum Gx: ', num2str(Gx_max)]);
disp(['Minimum Gx: ', num2str(Gx_min)]);
disp(['Maximum Gz: ', num2str(Gz_max)]);
disp(['Minimum Gz: ', num2str(Gz_min)]);
disp(['Angular velocity: ', num2str(w_max)]);

% Figure for g-force over time
figure;
plot(t, Gz, 'LineWidth', 2);
set(gca, 'FontSize', 15);
grid on;
xlabel('t', 'FontSize', 18);
ylabel('$G_z$', 'Interpreter', 'latex', 'FontSize', 18);
xlim([0, t3]);
legend('$G_z$', 'Interpreter', 'latex', 'FontSize', 18);

figure;
plot(t, Gx, 'LineWidth', 2);
set(gca, 'FontSize', 15);
grid on;
xlabel('t', 'FontSize', 18);
ylabel('$G_x$', 'Interpreter', 'latex', 'FontSize', 18);
xlim([0, t3]);
legend('$G_x$', 'Interpreter', 'latex', 'FontSize', 18);

% Perform FFT on Gx and Gz
N = length(t); % Number of samples
Fs = N / t3; % Sampling frequency (in Hz)

% Compute the magnitude of the FFT
Gx_fft = fft(Gx);
Gz_fft = fft(Gz);
Gx_magnitude = abs(Gx_fft / N);
Gz_magnitude = abs(Gz_fft / N);

% Compute the corresponding frequencies
f = (0:N-1) * (Fs / N); 

% Consider only the first half of the FFT (positive frequencies)
half_idx = 1:floor(N/2);
f_half = f(half_idx);
Gx_magnitude_half = Gx_magnitude(half_idx);
Gz_magnitude_half = Gz_magnitude(half_idx);

% Find the frequency with the maximum magnitude (excluding DC component)
[~, Gx_max_idx] = max(Gx_magnitude_half(2:end)); 
[~, Gz_max_idx] = max(Gz_magnitude_half(2:end));

% The most significant frequency for Gx and Gz
Gx_frequency = f_half(Gx_max_idx + 1); 
Gz_frequency = f_half(Gz_max_idx + 1);

% Display the results
disp(['Frequency of Gx: ', num2str(Gx_frequency)]);
disp(['Frequency of Gz: ', num2str(Gz_frequency)]);