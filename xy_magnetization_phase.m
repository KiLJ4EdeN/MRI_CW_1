B0 = 1.5e3; % 1.5 T = 1500 mT
% using units of mT and ms for the bloch_rotate function

% start at equilibrium
Mstart = [0,0,1].'; 
% after RF excitation
Mstart = [1,0,0].'; 


Bstatic = [0,0,B0];

% animate rotation?
dt = .1e-6; % .1 ns = .1e-6 ms
N = 300;
t = [1:N]*dt;
Mall = zeros(3,N);
Mall(:,1) = Mstart;
for It = 1:N-1
    Mall(:,It+1) = bloch_rotate(Mall(:,It),dt,Bstatic);
end
plot(t,Mall)
xlabel('time (ms)'), ylabel('Magnetization')
legend({'M_X', 'M_Y', 'M_Z'}, 'location', 'north'), legend boxoff
title(['Precession of the transverse magnetization'])


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