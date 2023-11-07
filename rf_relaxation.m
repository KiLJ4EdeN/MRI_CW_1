t = linspace(0,1); % s

M0=1;
T1 = .8;  T2 = .1; % s
M_equilibrium = [0,0,M0].';

flip= 90;

gammabar = 42.58; % kHz/mT
T = 1; % 1 ms pulse duration
% calculate RF pulse amplitude (milliTesla)

% first one for alpha_x, second one for alpha_y
B10 = flip*pi/180 / (2*pi*gammabar*T)
% B10 = (flip*pi/180 / (2*pi*gammabar*T)) * 1i

% apply RF tip
M_start = bloch_rftip(M_equilibrium, T, B10)

Mall = zeros(3,length(t));

for It = 1:length(t)
    Mall(:,It) = bloch_relax(M_start,t(It),M0,T1, T2);
end

plot(t,Mall)
xlabel('time (s)'), ylabel('Magnetization')
legend({'M_X', 'M_Y', 'M_Z'}, 'location', 'north'), legend boxoff
title(['Relaxation after a ' num2str(flip) '-degree flip angle'])


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