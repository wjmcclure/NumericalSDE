close all
k0 = 0.4; k1 = 0.01; k2 = 1; k3 = 1; k4 = 0.2; j3 = 0.05; j4 = 0.05;
S = @(t) 6*((t>20).*(t<30) + (t>50) + (t>60)).*(t<70);
G = @(u,v,J,K) 2*u*K*...
    (v-u+v*J+u*K + sqrt((v-u+v*J+u*K)^2 -4*(v-u)*u*K))^-1;
[T,U] = ode45( @(t,u) progressode(t,u,S,G,k0,k1,k2,k3,k4,j3,j4),[0 100],0);

figure('Position',[134 393 771 413]);
subplot(2,1,1);
plot(T,S(T),'Color',[0,0,1],'LineWidth',2);
axis([0 100 -1 13]);
title('One-way toggle response');
ylabel('S concentration');
subplot(2,1,2);
plot(T,U,'Color',[1,0,0],'LineWidth',2);
axis([0 100 -0.1, 0.6]);
xlabel('time');
ylabel('R concentration');

function dudt = progressode(t,u,S,G,k0,k1,k2,k3,k4,j3,j4)
    dudt = k0*G(k3*u,k4,j3,j4) + k1*S(t) - k2*u;
end
