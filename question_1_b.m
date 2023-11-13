% Constants
B0 = 1.5; % Tesla ( not used in rotating frame )
B1 = 0.05; % Gauss
B1 = B1 * 1e-04; % Gauss to Tesla
M0 = [0; 0; 1]; % Initial magnetization vector
B1_time = 7.35; % in ms

% Time range
tspan = [0, B1_time/1000]; % Time span in milliseconds

% Solve the ODE
[t, M] = ode45(@(t, M) bloch_equations(t, M, B0, B1), tspan, M0);

% initial and final values of magnetization
equilibrium_M = [M(1, 1), M(1, 2), M(1, 3)]
final_M = [M(end, 1), M(end, 2), M(end, 3)]

% Plot Mx, My, Mz
figure
plot(M(:, 1), 'r')
hold on
plot(M(:, 2), 'b')
plot(M(:, 3), 'm')
title("Magnetization values 2D")
legend('Mx','My', 'Mz')
pause(2)
saveas(gcf,'Figures/question_1_b_2d_magnetization.png')
close all

% Plot the magnetization trajectory all in one
figure
plot3(M(1, 1), M(1, 2), M(1, 3), 'ro', 'linewidth', 4)
hold on
plot3(M(end, 1), M(end, 2), M(end, 3), 'go', 'linewidth', 4)
plot3(M(:, 1), M(:, 2), M(:, 3), 'b');
xlabel('M_x');
ylabel('M_y');
zlabel('M_z');
title('Magnetization Vector Trajectory');
grid on;
pause(2)
saveas(gcf,'Figures/question_1_b_magnetization_trajectory.png')
close all

% figure,
view([1,1,1]);
am = 1; % default is 1
axis([-am am -am am -am am]);
hold on;
grid on;
f = gcf;
width = 800;
height = 800;
Pix_SS = get(0,'screensize');
set(f,'Position',[(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height])
set(f,'Renderer','ZBuffer')
% Quiver3 plot over time
writerObj = VideoWriter('Videos/question_1_b.avi');
open(writerObj);
for n=1:length(M)
    M_vector = M(n, :);
    Mx = M_vector(1);
    My = M_vector(2);
    Mz = M_vector(3);
    h1=quiver3(0,0,0,Mx,My,Mz, 'Color', 'r');
    h2=quiver3(0,0,0,Mx,My,0, 'Color', 'b');
    h3=quiver3(0,0,0,0,0,Mz, 'Color', 'm');
    F = getframe(f) ;           %// Capture the frame
    writeVideo(writerObj,F)  %// add the frame to the movie
    pause(0.1)
    delete(h1)
    delete(h2)
    delete(h3)
end
close all
close(writerObj);

% Bloch equations
function dMdt = bloch_equations(t, M, B0, B1)
    % Define stuff again here, globals cant be used
    off_resonance = 1e3; % off resonance
    gamma = 42.6e6; % Gyromagnetic ratio for protons in Hz/T
    B1 = [B1; 0; 0]; % converted B1 to rotating frame from book (B1i)
    w0 = (gamma + off_resonance) * B0;
    w_rf = gamma * B0;
    delta_w = abs(w0 - w_rf);
    B_eff = B1 + (B0 - delta_w) / gamma * [0; 0; 1]; % because of off resonance eq 
    % eq 3.75 book
    dMdt = gamma * cross(M, B_eff); % Bloch equations in rotating frame
end
