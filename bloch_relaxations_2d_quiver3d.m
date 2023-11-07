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

[Rows, Columns] = size(Mall);

figure,
view([1,1,1]);
am = 1; % default is 1
axis([-am am -am am -am am]);
hold on;
grid on;

for k=1:Columns
    Mx = Mall(1, k);
    My = Mall(2, k);
    Mz = Mall(3, k);
    h1=quiver3(0,0,0,Mx,My,Mz, 'Color', 'r');
    h2=quiver3(0,0,0,Mx,My,0, 'Color', 'b');
    h3=quiver3(0,0,0,0,0,Mz, 'Color', 'm');
    pause(0.3)
    delete(h1)
    delete(h2)
    delete(h3)
end

close all

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