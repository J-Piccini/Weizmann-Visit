close all;
clear;
clc;
%% Environment parameter

tmax = 2e3;    % Simulation length
dt = 1;    % Time step
r = 0.1;    % Influence radius

Ly = 2;    % Height of the environment
Lx = Ly * 7;    % Lx / Ly is taken from experimentla set up (2mm thickness and 1.4cm diameter)
Lgridy = 1 : r : Ly;    % y-grid (for speed purposes)
Lgridx = 1 : r : Lx;    % x-grid (for speed purposes)
%% Species characterisation

noise0 = 0.1;    % Standard value for noise
knoise = 1;    % Noise modifier (for second species)
noise1 = noise0;    % Noise for species 1 
noise2 = knoise * noise0;    % Noise for species 2 

N = 300;    % Number of particles
N1 = round(0.5 * N);    % 50-50 particles division, round for cases when N is odd
N2 = N - N1;

vel = 0.03;    % Standard value for velocity
kvel = 2;    % Velocity modifier parameter (for second species)
velv = [vel * ones(1, N1), kvel * vel * ones(1, N2)];    % Velocity vector
%%

x = Lx * rand(1, N);
y = Ly * rand(1, N);
theta = 2 * pi * (rand(1, N) - 0.5);

MX = zeros(N, tmax);
MY = zeros(N, tmax);
MTH = zeros(N, tmax);

for time = 1 : tmax

    % Utilies generation
    segloc_x = zeros(1, N);
    segloc_y = zeros(1, N);
    Mnear = cell(1, N);
    MM = zeros(N);

    % Boundaries check
    x(x < 0) = -x(x < 0);
    x(x > Lx) = 2 * Lx - x(x > Lx);
    y(y < 0) = -y(y < 0);
    y(y > Ly) = 2 * Ly - y(y > Ly);

    % Box of each particle
    for i = 1 : N 
        idx_x = x(i) > Lgridx;
        segloc_x(i) = min(find(idx_x == 0));

        idx_y = y(i) > Lgridy;
        segloc_y(i) = min(find(idx_y == 0));
    end


    % Particles in the interaction radius
    for i = 1 : N 
        idx_i = segloc_x(i);
        idy_i = segloc_y(i);

        %             cand_x = find(segloc_x(i + 1 : end) >= idx_i - 1 & segloc_x(i + 1 : end) <= idx_i + 1) + i;
        %             cand_y = find(segloc_y(i + 1 : end) >= idy_i - 1 & segloc_y(i + 1 : end) <= idy_i + 1) + i;
        %             This is more efficient, however, there is a need to update
        %             near_xy = intersect(cand_x, cand_y);
        %             Mnear{i} = near_xy;

        control_x = find(segloc_x(1 : end) >= idx_i - 1 & segloc_x(1 : end) <= idx_i + 1);
        control_y = find(segloc_y(1 : end) >= idy_i - 1 & segloc_y(1 : end) <= idy_i + 1);
        near2 = intersect(control_x, control_y);

        distances = sqrt((x(i) - x(near2)).^2 + (y(i) - y(near2)).^2);

        Mnear{i} = near2(distances <= r);

        if ~isempty(Mnear{i})
            avg_th(i) = atan2(mean(sin(theta(Mnear{i}))), mean(cos(theta(Mnear{i}))));
        else
            avg_th(i) = theta(i);
        end
    end

    %% Update position & Boundaries Implementation

    % x(t + 1) = x(t) + v(t) * delta_t
    x0 = x;
    y0 = y;
    % First update to find which particle bounce
    x = x + velv .* cos(theta) * dt;
    y = y + velv .* sin(theta) * dt;  

    noise_v = [noise1 * (rand(1, N1) - 0.5), noise2 * (rand(1, N2) - 0.5)];        
    theta = avg_th + noise_v;

    % Bouncing angles
    theta(x < 0 | x > Lx) = pi - theta(x < 0 | x > Lx); 
    theta(y < 0 | y > Ly) = - theta(y < 0 | y > Ly);

    % True update
    x = x0 + velv .* cos(theta) * dt;
    y = y0 + velv .* sin(theta) * dt;  

    % Trajectories storing
    for i = 1 : N
        MX(i, time) = x(i);
        MY(i, time) = y(i);
    end
    % Orientation storing
    MTH(:, time) = theta(:);

    %% Plotting at each time interval 

%     plot(x(1 : N1), y(1 : N1), '.', x(N1 + 1 : N), y(N1 + 1 : N), '.', 'MarkerSize', 15)
%     xlim([0, Lx])
%     ylim([0, Ly])
% 
%     pause(vel)
end
%% Final distribution plot

figure;
hold on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 24)
plot(x(1 : N1), y(1 : N1), '.', x(N1 + 1 : N), y(N1 + 1 : N), '.', 'MarkerSize', 15)
xlim([0, Lx])
ylim([0, Ly])
%% Further plots

additional = true;

if additional

    id1 = ceil(N1 * rand);
    id2 = N1 + ceil(N2 * rand);
    figure;
    hold on
    grid on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 24)

    plot(MX(id1, :), MY(id1, :), 'LineWidth', 1.1)
    plot(MX(id2, :), MY(id2, :), 'LineWidth', 1.1)

    plot(MX(id1, 1), MY(id1, 1), 'p', 'LineWidth', 1.1)
    plot(MX(id1, end), MY(id1, end), 'o', 'LineWidth', 1.1)

    plot(MX(id2, 1), MY(id2, 1), 'p', 'LineWidth', 1.1)
    plot(MX(id2, end), MY(id2, end), 'o', 'LineWidth', 1.1)

    legend('Species 1', 'Species 2', 'Species 1 -Starting point', 'Species 1 - Ending point', ...
        'Species 2 - Starting point', 'Species 2 - Ending point', ...
        'Interpreter', 'latex', 'FontSize', 24, 'Location', 'best')
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24)
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 24)

    axis tight
end