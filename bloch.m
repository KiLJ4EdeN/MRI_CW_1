function bloch(T1,T2,f)
t=linspace(0,3*T1,10000);
Mz=(1-exp(-t/T1));
Mx=exp(-t/T2).*(cos(2*pi.*f.*t)+sin(2*pi*f.*t));
My=exp(-t/T2).*(cos(2*pi*f.*t)-sin(2*pi*f.*t)); 
figure,
view([1,1,1]);
axis([-max(Mx) max(Mx) -max(My) max(My) 0 max(Mz)]);
grid on
hold on
for i=1:length(t)
% plot3(Mx(i),My(i),Mz(i),'dg','MarkerSize',5);
% h2=plot3([0,Mx(i)],[0,My(i)],[0,Mz(i)],'r','LineWidth',2);
% h3=plot3([0,Mx(i)],[0,My(i)],[0,0],'b','LineWidth',2);
% h4=plot3([0,0],[0,0],[0,Mz(i)],'m','LineWidth',2);
% quiver3(0,0,0,Mx(i),My(i),Mz(i),'dg','MarkerSize',5);
h2=quiver3(0,0,0,Mx(i),My(i),Mz(i), 'Color', 'r');
h3=quiver3(0,0,0,Mx(i),My(i),0, 'Color', 'b');
h4=quiver3(0,0,0,0,0,Mz(i), 'Color', 'm');
title(strcat(num2str(t(i)),'[s]'))
drawnow limitrate
delete(h2)
delete(h3)
delete(h4)
end
end