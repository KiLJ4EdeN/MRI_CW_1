t = linspace(0,1); % s

T1 = .8;  T2 = .1; % s
M0 = 1;

% M_start = [M0, 0, 0].'; % start along MX
% M_start = [0, M0, 0].'; % start along MY
M_start = [M0/3, M0/2, 0].'; % start along MXY

Mall = zeros(3,length(t));

for It = 1:length(t)
    Mall(:,It) = bloch_relax(M_start,t(It),M0,T1, T2);
end

plot(t,Mall)
xlabel('time (s)'), ylabel('Magnetization')
legend({'M_X', 'M_Y', 'M_Z'}, 'location', 'north'), legend boxoff
title(['Relaxation'])



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