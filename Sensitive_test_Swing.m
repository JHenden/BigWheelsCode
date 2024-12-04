%% chebivp.m -- Solving for different combinations of v and d
% Automatically modified to solve for specific values of r, M, v, and d.

%% Problem description.
% Solving
%   u'' + sin(u) - cos(u)*sin(u+4/pi) = 0,
% for t in [0 10], subject to
%   u'(0) = 0,
%   u(0) = 0.

% Define parameter values
v_values = [10, 15, 20, 25, 30, 35];
d_values = [0, 45, 90, 135, 180, 225, 270, 315]; % d_values repersent the set of degrees

% Define constants
M = 200;
l = 1.5;
A0 = 4 * l^2;
g = 9.81;
Cd = 1.225;
ro = 1.29;
t1 = 0.5;
t2 = 1.5;
t3 = 3;

% Define the domain
dom = [0 10];

% Function to define vw as a chebfun
function vw = get_vw(v_max, t1, t2, t3, dom)
    vw = chebfun({...
        @(t) (v_max / t1) * t, ...             % From 0 to t1: Linear increase
        @(t) v_max, ...                         % From t1 to t2: Stay v_max
        @(t) v_max * (1 - (t - t2) / (t3 - t2)), ... % From t2 to t3: Linear decrease
        @(t) 0}, ...                           % From t3 to dom(2): Zero velocity
        [0, t1, t2, t3, dom(2)]);              % Domain breakpoints for each piece
end

% Loop over all combinations of v and d
for i = 1:length(v_values)
    for j = 1:length(d_values)
        v_max = v_values(i);
        d = d_values(j);
        angle = d * pi / 180; % d repersent the degree, where rad = d * pi / 180
        
        % Define vw as a chebfun using function
        vw = get_vw(v_max, t1, t2, t3, dom);
        
        % Assign the differential equation to a chebop on that domain.
        N = chebop(@(t, u) M * l * diff(u, 2) + M * g * sin(u) - 0.5 * sqrt(2) * A0 * Cd * ro * (vw.^2) .* cos(u-angle) .* sin((u + 4 ./ pi)), dom);
        
        % Set up the rhs of the differential equation so that N(u) = rhs.
        rhs = 0;
        
        % Assign boundary conditions to the chebop.
        N.lbc = @(u) [diff(u); u]; % Corrected to explicitly set boundary conditions.
        
        % Setup preferences for solving the problem.
        options = cheboppref();
        options.ivpSolver = 'ode113';
        
        %% Solve!
        try
            u = solveivp(N, rhs, options);
            
            % Calculate w and its derivative
            w = diff(u);
            diff_w = diff(w);
            
            % Calculate gforce using element-wise operations
            ar = l * (w).^2;
            at = l * diff_w;
            Gx = sin(u) + at / g;
            Gz = cos(u) + ar / g;
            
            % Calculate the maximum and minimum g-force
            Gx_max = max(Gx);
            Gx_min = min(Gx);
            Gz_max = max(Gz);
            Gz_min = min(Gz);
            
            % Calculate the maximum velocity
            v = w .* l;
            vmax = max(v);
            u_max = max(u);
            
            % Display the results for this combination (vertically)
            fprintf('Results for v_max = %d, d = %d degrees\n', v_max, d);
            fprintf('  Maximum Gx: %.4f\n', Gx_max);
            fprintf('  Minimum Gx: %.4f\n', Gx_min);
            fprintf('  Maximum Gz: %.4f\n', Gz_max);
            fprintf('  Minimum Gz: %.4f\n', Gz_min);
            fprintf('  Maximum velocity: %.4f\n', vmax);
            disp(['  beta: ', num2str(u_max)]);
            fprintf('----------------------------------\n');
        catch ME
            % If an error occurs during solving, display a warning
            fprintf('Failed to solve for v_max = %d, d = %d degrees: %s\n', v_max, d, ME.message);
        end
    end
end
