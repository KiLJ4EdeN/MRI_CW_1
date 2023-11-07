% lab frame

gammabar = 42.58; % kHz/mT

B0 = 10; %  10 mT for simplicity to visualize rotation
f0 = gammabar*B0 % kHz

M0 = 1;
M_equilibrium = [0,0,M0].';

% RF pulse paarameters
T_RF =  1; % ms
t = linspace(0, T_RF, 4000);

RF_flip_angle = pi/2; % radians
B10 = RF_flip_angle / (2*pi*gammabar*T_RF) % mT

% RF not applied at resonance frequency (constant magnetic field in X)
% lab frame
B = [B10;0;B0];  

Mall = zeros(3,length(t));
Mall(:,1) = M_equilibrium;
for It = 1:length(t)-1
    Mall(:,It+1) = bloch_rotate(Mall(:,It),t(It+1) - t(It),B);
end

% RF not applied at resonance frequency (constant magnetic field in X)
% lab frame
B = [B10;0;B0];  

Mall = zeros(3,length(t));
Mall(:,1) = M_equilibrium;
for It = 1:length(t)-1
    Mall(:,It+1) = bloch_rotate(Mall(:,It),t(It+1) - t(It),B);
end

plot(t,Mall)
xlabel('time (ms)'), ylabel('Magnetization')
legend({'M_X', 'M_Y', 'M_Z'}, 'location', 'north'), legend boxoff
title(['Non-resonant magnetic field (in RF direction)'])



% RF pulse at Larmor frequency
% lab frame
B = [B10*cos(2*pi*f0*t);B10*-sin(2*pi*f0*t);B0*ones(1,length(t))];  

Mall = zeros(3,length(t));
Mall(:,1) = M_equilibrium;
for It = 1:length(t)-1
    Mall(:,It+1) = bloch_rotate(Mall(:,It),t(It+1) - t(It),B(:,It));
end
figure
plot(t,Mall)
xlabel('time (ms)'), ylabel('Magnetization')
legend({'M_X', 'M_Y', 'M_Z'}, 'location', 'northwest'), legend boxoff
title(['Resonant RF pulse'])


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
