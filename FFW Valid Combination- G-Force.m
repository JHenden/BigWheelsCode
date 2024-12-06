%% Valid T and r combinations for Fast Ferris Wheels

% Define parameters
w_values = linspace(0.70, 1.25, 1000);  % Possible values for period time T
r_values = linspace(5, 30, 1000);  % Possible values for radius r
g = 9.81;                       % Gravitational acceleration in m/s^2

% Define the parameters of the piecewise function
t1 = 25;
t2 = 50;
t3 = 75;
t = linspace(0, t3, 1000);  

% Initialize an array to store valid T and r combinations
valid_combinations = [];
labels = [];

% Loop over all combinations of T and r
for w_max = w_values
    for r = r_values
        % Calculate maximum angular velocity
        
        % Define the piecewise function for angular velocity over time
        f = @(t) (t >= 0 & t < t1) .* (w_max / t1 .* t) + ...
                 (t >= t1 & t < t2) .* w_max + ...
                 (t >= t2 & t <= t3) .* (w_max - (w_max / (t3 - t2)) .* (t - t2));
        
        % Define the piecewise function for derivative angular velocity over time
        h = @(t) (t >= 0 & t < t1) .* (w_max / t1) + ...
                 (t >= t1 & t < t2) .* 0 + ...
                 (t >= t2 & t <= t3) .* (- (w_max / (t3 - t2)));
        % Calculate the angular velocity and its derivative at each time point
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
        
        % Check if exactly one g-force value meets the specified limit condition, while the others are below their respective limits
        if (abs(Gx_max - 6) < 0.001 && Gx_min > -2 && Gz_max < 6 && Gz_min > -2) || ...
           (abs(Gx_min + 2) < 0.001 && Gx_max < 6 && Gz_max < 6 && Gz_min > -2) || ...
           (abs(Gz_max - 6) < 0.001 && Gx_max < 6 && Gx_min > -2 && Gz_min > -2) || ...
           (abs(Gz_min + 2) < 0.001 && Gx_max < 6 && Gx_min > -2 && Gz_max < 6)
            % Store the valid combination of T and r
            valid_combinations = [valid_combinations; w_max, r];
        end
    end
end

% Display the valid combinations of T and r
if isempty(valid_combinations)
    disp('No valid combinations of w_max and r found.');
else
    disp('Valid combinations of w_max and r:');
    disp(array2table(valid_combinations, 'VariableNames', {'w_max', 'r'}));
    
    % Plot the valid combinations
    figure;
    scatter(valid_combinations(:, 1), valid_combinations(:, 2), 20, 'filled');
    xlabel('Angular velocity $\omega_{\theta}$ (rad/s)', 'Interpreter', 'latex', 'FontSize', 28);
    ylabel('Radius r (m)', 'FontSize', 28);
    set(gca, 'FontSize', 18);
    grid on;
end
