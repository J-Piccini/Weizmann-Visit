close all;
clear;
clc;
%% Environment definition

N = 100;
density = 3;
Ly = floor(sqrt(N / (density * 4)));
Lx = 4 * Ly;
Lx = 32;
Ly = 32;
N = Lx * Ly * density;

r0 = 1;
grid_x = linspace(0, Lx, ceil(Lx / r0));
grid_y = linspace(0, Ly, ceil(Ly / r0));

tSim = 1e3;
dt = 1;
%%

x = Lx * rand(1, N);
y = Ly * rand(1, N);
theta = 2 * pi * (rand(1, N) - 0.5);

vel0 = 0.5;
noise = 0.3;
loop = true;
time = 0;
while loop
    time = time + 1;
    

    x(x <= 0) = Lx + x(x <= 0);
    x(x >= Lx) = x(x >= Lx) - Lx;

    y(y >= Ly) =  2 * Ly - y(y >= Ly);
    y(y <= 0) = abs(y(y <= 0));

    % Utilies generation
    segloc_x = zeros(1, N);
    segloc_y = zeros(1, N);

    % Particles in the interaction radius
    for i = 1 : N
        idx_x = x(i) > grid_x;
        segloc_x(i) = min(find(idx_x == 0)) - 1;

        idx_y = y(i) > grid_y;
        segloc_y(i) = min(find(idx_y == 0)) - 1;
    end

    clearvars idx_x idx_y
    
    for i = 1 : N
        tX = segloc_x(i);
        nX = find(segloc_x >= tX - 1 & segloc_x <= tX + 1);
        if nX
            tY = segloc_y(i);
            nY = find(segloc_y >= tY - 1 & segloc_y <= tY + 1);

            if nY
                cand = intersect(nX, nY);

                distances = sqrt((x(i) - x(cand)).^2 + (y(i) - y(cand)).^2);
                near = cand(distances <= r0);

                if near
                    avg_th(i) = atan2(mean(sin(theta(near))), mean(cos(theta(near))));
                else
                    avg_th(i) = theta(i);
                end
            end
        end
    end

    theta = avg_th + 2 * noise * pi * (rand(1, N) - 0.5);
    x = x + vel0 * cos(theta) * dt;
    y = y + vel0 * sin(theta) * dt;

    %% Reflective Boundary

%     theta(x < 0 | x >= Lx) = pi - theta(x < 0 | x >= Lx);
%     theta(y < 0 | y >= Ly) = - theta(y < 0 | y >= Ly);
% 
%     x(x >= Lx) =  2 * Lx - x(x >= Lx);
%     x(x <= 0) = abs(x(x <= 0));

%     y(y >= Ly) =  2 * Ly - y(y >= Ly); % Reflective y
%     y(y <= 0) = abs(y(y <= 0));

    %% Periodic Boundary
    x(x < 0) = Lx + x(x < 0); % Reflective x
    x(x > Lx) = x(x > Lx) - Lx;

    y(y < 0) = Ly + y(y < 0);
    y(y > Ly) = y(y > Ly) - Ly;

    
    %toc
    
    if ~mod(time, 10)
        fprintf('\nTime: %d', time)
        if ~mod(time, 100)
            clc;
        end
        
    end
    plot(x, y, '.', 'MarkerSize', 15)
    hold on
    quiver(x, y, vel0 * cos(theta), vel0 * sin(theta), 0);
    xlim([0 Lx]);
    ylim([0 Ly]);
    hold off
    pause(0.001)
end
