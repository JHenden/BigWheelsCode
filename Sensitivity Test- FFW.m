%% Sensitive test for Fast Ferris Wheels

% Define parameters
T_values = [5,6,7,8,9];                  % Period times
r_values = [5,8,11,14,17,20];           % Radius values
g = 9.81;                                   % Gravitational acceleration in m/s^2

% Define the parameters of the piecewise function
t1 = 25;                         
t2 = 50;                        
t3 = 75;                        

% Define the time range for calculating and plotting
t = linspace(0, t3, 1000);      % Time vector

% Loop through all combinations of T and r
for i = 1:length(T_values)
    for j = 1:length(r_values)
        % Current values of T and r
        T = T_values(i);
        r = r_values(j);
        
        % Calculate maximum angular velocity
        w_max = 2 * pi / T;
        v_max = r * w_max;
        
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

        % Display the maximum and minimum g-force in the command window
        disp(['For T = ', num2str(T), ', r = ', num2str(r), ':']);
        disp(['Maximum Gx: ', num2str(Gx_max)]);
        disp(['Minimum Gx: ', num2str(Gx_min)]);
        disp(['Maximum Gz: ', num2str(Gz_max)]);
        disp(['Minimum Gz: ', num2str(Gz_min)]);
        disp(['Maximum v: ', num2str(v_max)]);
        disp('------------------------------------');
    end
end
