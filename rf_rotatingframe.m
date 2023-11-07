% lab frame

gammabar = 42.58; % kHz/mT

B0 = 10; %  10 mT for simplicity to visualize rotation
f0 = gammabar*B0; % kHz

M0 = 1;
M_equilibrium = [0,0,M0].';

% RF pulse paarameters
T_RF =  1; % ms
t = linspace(0, T_RF, 4000);

RF_flip_angle = pi/2; % radians
B10 = RF_flip_angle / (2*pi*gammabar*T_RF); % mT

% rotating frame
% RF pulse at Larmor frequency
B = [B10;0;0];  

Mall = zeros(3,length(t));
Mall(:,1) = M_equilibrium;
for It = 1:length(t)-1
    Mall(:,It+1) = bloch_rotate(Mall(:,It),t(It+1) - t(It),B);
end

plot(t,Mall)
xlabel('time (ms)'), ylabel('Magnetization')
legend({'M_X', 'M_Y', 'M_Z'}, 'location', 'north'), legend boxoff
title(['Resonant RF pulse in rotating frame'])

