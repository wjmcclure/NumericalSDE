close all
a = 0.7; %Stability parameter

n = 1e4;
Tfinal = 100;
T = linspace(0,Tfinal,n);
h = Tfinal / n;

%Sniffer parameters
k1 = 50;
k2 = 13;
k3 = 1;
k4 = 1;
%Goldbeter-Koshland function
G = @(u,v,J,K) 2*u*K/(v-u+v*J+u*K + sqrt((v-u+v*J+u*K)^2-4*(v-u)*u*K));
%Breaker/fuse parameters
K0 = 0.4;
K1 = 0.01;
K2 = 1;
K3 = 1;
K4 = 0.2;
J3 = 0.05;
J4 = 0.05;

U = zeros(4,n); %U = (S,R,X,Rb)'
Sbase = 0.5; %Baseline signal level
Rss = k1*k4 / (k2*k3); %Steady-state value of R
U(:,1) = [ Sbase; Rss; 0; 0.04];

%It appears that simulation with negative signal is unstable.
rng(45350987);
Unif = rand(n,1);
P = makedist('Stable','alpha', a, 'beta', 0, 'gam', 1, 'delta', 0);
SaS = icdf(P,Unif);
P = makedist('normal','mu',0,'sigma',1);
Normal = icdf(P,Unif);
%fprintf('Highest Levy noise is %.2e\n',max(SaS));
%eps = h^a*(3e3*rand +1) / (2*max(SaS)); %limit jumps in S
eps = 1e-2;
delta = 4e-2;
idx = find(T>=20,1);
times = [10 10 55]; %Nothing -- linear -- brownian -- levy
for i=2:n
%     if mod(i,50)==0
%         fprintf('%d\n',i);
%     end
    if T(i) < times(1)
        U(1,i) = U(1,i-1);
    else if T(i) < times(2)
        U(1,i) = U(1,i-1) + h*8e-2;
    else if T(i) < times(3) 
        U(1,i) = U(1,i-1) + h^(1/2)*delta*Normal(i);
        lasti = i;
    else
        U(1,i) = U(1,i-1) + eps*h^(1/a)*SaS(i);
    end
    end
    end
    U(2,i) = U(2,i-1) + h*( k1*U(1,i-1) - k2*U(2,i-1)*U(3,i-1) );
    U(3,i) = U(3,i-1) + h*( k3*U(1,i-1) - k4*U(3,i-1) );
    U(4,i) = U(4,i-1) + K0*G(K3*U(4,i-1),K4,J3,J4) + K1*U(2,i-1) - ...
                K2*U(4,i-1);
end

figure('Position',[134 393 771 620]); 
subplot(3,1,1);
plot(T,U(1,:),'Color',[0,0,1],'LineWidth',2);
ylabel('S concentration');
title('Sniffer + one-way with Brownian and Levy noise');
dist = abs( min(U(1,:)) - max(U(1,:)) );
axis([times(1) Tfinal min(U(1,:))-0.1*dist max(U(1,:))+0.1*dist]);
subplot(3,1,2);
plot(T,U(2,:),'Color',[1,0,0],'LineWidth',2);
ylabel('R1 concentration');
dist = abs( min(U(2,:)) - max(U(2,:)) );
axis([times(1) Tfinal min(U(2,:))-0.1*dist max(U(2,:))+0.1*dist]);
subplot(3,1,3);
plot(T,U(4,:),'Color',[0,0.4,0],'LineWidth',2);
xlabel('time');
ylabel('R2 concentration');
dist = abs( min(U(4,:)) - max(U(4,:)) );
axis([times(1) Tfinal min(U(4,:))-0.1*dist max(U(4,:))+0.1*dist]);