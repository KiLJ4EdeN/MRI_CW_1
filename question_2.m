% A) lookup: Rx.m, Ry.m, Rz.m
% B) lookup the bloch_relax function in this file
% C)
% divide into thousand so each t is a ms
t = linspace(0, 1, 1000); % s
% total Magnetization at eq
M0 = 1;
% define T1 T2 in s
T1 = 1.;  T2 = .1; % s
% define M
M_equilibrium = [0, 0, M0].';

gammabar = 42.58; % kHz/mT
T = 1; % 1 ms pulse duration

% calculate RF pulse amplitude (milliTesla)
% desired flip angle
flip = 60;

% first one for alpha_x, second one for alpha_y
% calculate B1 based on the formula
B10 = flip*pi/180 / (2*pi*gammabar*T);
% B10 = (flip*pi/180 / (2*pi*gammabar*T)) * 1i

Magnetization_history = {};
Magnetization_cmt = {};
Magnetization_history{1} = M_equilibrium;
Magnetization_cmt{1} = "M equilibrium";

% apply RF tip
M_after_60_x = bloch_rftip(M_equilibrium, T, B10)
Magnetization_history{2} = M_after_60_x;
Magnetization_cmt{2} = "M after 60 x";

relaxation_time_1 = 15; %ms

M_rel_1 = zeros(3, relaxation_time_1);

for It = 1:relaxation_time_1
    M_rel_1(:,It) = bloch_relax(M_after_60_x,t(It),M0,T1, T2);
end

Magnetization_history{3} = M_rel_1;
Magnetization_cmt{3} = "M during relaxation 1";

% see the last magnetization vector
M_after_relaxation_1 = M_rel_1(:, end)

Magnetization_history{4} = M_after_relaxation_1;
Magnetization_cmt{4} = "M after relaxation 1";

% B1 number 2
flip = 45;
B10 = flip*pi/180 / (2*pi*gammabar*T) + 1i;

M_after_45y = bloch_rftip(M_equilibrium, T, B10)

Magnetization_history{5} = M_after_45y;
Magnetization_cmt{5} = "M after 45y";

relaxation_time_2 = 10; %ms

M_rel_2 = zeros(3, relaxation_time_2);

for It = 1:relaxation_time_2
    M_rel_2(:,It) = bloch_relax(M_after_45y,t(It),M0,T1, T2);
end

Magnetization_history{6} = M_rel_2;
Magnetization_cmt{6} = "M during relaxation 2";

% see the last magnetization vector
M_after_relaxation_2 = M_rel_2(:, end)

Magnetization_history{7} = M_after_relaxation_2;
Magnetization_cmt{7} = "M after relaxation 2";

figure,
view([1,1,1]);
am = 1; % default is 1
axis([-am am -am am -am am]);
hold on;
grid on;

for k=1:length(Magnetization_history)
    hist = Magnetization_history{k};
    [Rows, Columns] = size(hist);
    cmt = Magnetization_cmt{k};
    % is singular
    if (Columns == 1)
        Mx = hist(1);
        My = hist(2);
        Mz = hist(3);
        h1=quiver3(0,0,0,Mx,My,Mz, 'Color', 'r');
        h2=quiver3(0,0,0,Mx,My,0, 'Color', 'b');
        h3=quiver3(0,0,0,0,0,Mz, 'Color', 'm');
        title(cmt)
        pause(2)
        delete(h1)
        delete(h2)
        delete(h3)
    end
end

close all;

figure,
view([1,1,1]);
am = 1; % default is 1
axis([-am am -am am -am am]);
hold on;
grid on;

for k=1:length(Magnetization_history)
    hist = Magnetization_history{k};
    [Rows, Columns] = size(hist);
    cmt = Magnetization_cmt{k};
    if (Columns ~= 1)
        for n=1:Columns
            M_vector = hist(:, n);
            Mx = M_vector(1);
            My = M_vector(2);
            Mz = M_vector(3);
            h1=quiver3(0,0,0,Mx,My,Mz, 'Color', 'r');
            h2=quiver3(0,0,0,Mx,My,0, 'Color', 'b');
            h3=quiver3(0,0,0,0,0,Mz, 'Color', 'm');
            title(cmt)
            pause(1)
            delete(h1)
            delete(h2)
            delete(h3)
        end
        pause(1)
    end
end

close all;

function [Mend] = bloch_rotate(Mstart, T, B)
% bloch_rotate - compute the rotation of the net magnetization for a given magnetic field
%
% INPUTS
%	Mstart - initial magnetization
%	T - duration [ms]
%	B = [Bx, By, Bz] - magnetic field [mT]
% OUTPUTS
%   Mend - final magnetization
%
% T and B can also be in units of [us] and [T], respectively

GAMMA = 42.58; % kHz/mT

flip = 2*pi*GAMMA * norm(B) * T;

eta = acos(B(3) / (norm(B)+eps));

theta = atan2(B(2), B(1));

Mend = Rz(-theta)*Ry(-eta)*Rz(flip)*Ry(eta)*Rz(theta)* Mstart;

end

function [Mend] = bloch_relax(Mstart, T, M0, T1, T2)
% bloch_relax - compute the effect of relaxation on the net magnetization
%
% INPUTS
%	Mstart - initial magnetization
%	T - duration [ms]
%	M0 - equilibrium magnetization (default = 1)
%	T1 - longitudinal relaxation time [ms]
%	T2 - transverse relaxation time [ms]
% OUTPUTS
%   Mend - final magnetization

Arelax = [exp(-T/T2) 0 0; ...
          0 exp(-T/T2) 0; ...
          0 0 exp(-T/T1)];
brecover = [0; 0; M0*(1-exp(-T/T1))];
	
Mend = Arelax*Mstart + brecover;
end

function [Mend] = bloch_rftip(Mstart, T, B1)
% bloch_rftip - compute the rotation due to RF (B1) on the net magnetization
%   in the rotating frame (neglecting effects of B0 and demodulating at the Larmor frequency)
%
% INPUTS
%	Mstart - initial magnetization
%	T - duration [ms]
%	B1 - RF amplitude, B1X+iB1Y [mT]
% OUTPUTS
%   Mend - final magnetization


Mend = bloch_rotate(Mstart, T, [real(B1) imag(B1), 0]);
end