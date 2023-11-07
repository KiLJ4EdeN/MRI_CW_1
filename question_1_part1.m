% Constants
B1 = 0.05;  % Gauss
B0 = 1.5;  % Tesla
M0 = [0; 0; 1];  % Initial magnetization vector

% Parameters
gamma_val = 2*pi*42.58e6;  % Gyromagnetic ratio (Hz/T)
omega = gamma_val*B0;  % Larmor frequency (Hz)

% Define the ODE system
dydt = @(t, M) gamma_val*B0*cross(M, [0; 0; 1]) + gamma_val*B1*cos(omega*t)*[1; 0; 0] - gamma_val*B1*sin(omega*t)*[0; 1; 0];

% Solve the ODE system
[t, M] = ode45(dydt, [0, 2*pi/omega], M0);

% Display the rotating coordinate system
figure;
plot3(M(:, 1), M(:, 2), M(:, 3));
xlabel('M_x');
ylabel('M_y');
zlabel('M_z');
title('Magnetization Vector Trajectory');
grid on;