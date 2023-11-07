t = linspace(0, 1, 10000); % s

M0 = 1;
T1 = .8;  T2 = .1; % s
M_equilibrium = [0, 0, M0].';

flip = 90;

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
    Mall(:,It) = bloch_relax_modified(M_start,t(It),M0,T1, T2, gammabar);
end

Mx = Mall(1, :);
My = Mall(2, :);
Mz = Mall(3, :);

% figure,
% view([1,1,1]);
% axis([-max(Mx) max(Mx) -max(My) max(My) 0 max(Mz)]);
% grid on
% hold on
% for i=1:length(t)
% plot3(Mx(i),My(i),Mz(i),'dg','MarkerSize',5);
% h2=plot3([0,Mx(i)],[0,My(i)],[0,Mz(i)],'r','LineWidth',2);
% h3=plot3([0,Mx(i)],[0,My(i)],[0,0],'b','LineWidth',2);
% h4=plot3([0,0],[0,0],[0,Mz(i)],'m','LineWidth',2);
% title(strcat(num2str(t(i)),'[s]'))
% drawnow limitrate
% delete(h2)
% delete(h3)
% delete(h4)
% end


quiver3(t,t,t,Mx,My,Mz)


function [Mend] = bloch_relax_modified(Mstart, T, M0, T1, T2, f)
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

Arelax = [exp(-T/T2).*(cos(2*pi.*f.*T)+sin(2*pi*f.*T)) 0 0; ...
          0 exp(-T/T2).*(cos(2*pi*f.*T)-sin(2*pi*f.*T)) 0; ...
          0 0 exp(-T/T1)];
brecover = [0; 0; M0*(1-exp(-T/T1))];
	
Mend = Arelax*Mstart + brecover;
% Mz = (1-exp(-t/T1));
% Mx = exp(-t/T2).*(cos(2*pi.*f.*t)+sin(2*pi*f.*t));
% My = exp(-t/T2).*(cos(2*pi*f.*t)-sin(2*pi*f.*t)); 
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