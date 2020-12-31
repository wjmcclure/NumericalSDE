%integration time and number of iterations
Tfin = 20;
n = 1e6;

%model parameters
s = 10;
r = 28;
b = 8/3;
eps = 0.05;

%drift and diffusion
b = @(X) [ -s*X(1)+s*X(2); r*X(1)-X(2)-X(1)*X(3); -b*X(3)+X(1)*X(2) ];
sig = @(X) sqrt(eps)*diag([X(1),X(2),X(3)]);

%setup
h = Tfin / n;
U = zeros(3,n); 
U(:,1) = [7.6; 6.2; 30.5];
T = linspace(0,Tfin,n);

%Euler-Maruyama method
for i=1:n-1
    U(:,i+1) = U(:,i) + b(U(:,i))*h + sig(U(:,i))*sqrt(h)*randn(3,1);
end

% %Normal boring picture
% plot3( U(1,:),U(2,:),U(3,:));
% view([40 20]);

% Pretty picture that takes a few seconds
Nc = 1e4;
f = figure;
for i=1:Nc
    v = 1+(i-1)*100:i*100;
    plot3( U(1,v),U(2,v),U(3,v),'Color',hsv2rgb([i/Nc,0.6,1]) );
    hold on
end
f.Color = 0.2*ones(1,3);
a = f.CurrentAxes;
a.Color = f.Color;
a.XColor = f.Color;
a.YColor = f.Color;
a.ZColor = f.Color;
a.XTick = [];
a.YTick = [];
a.ZTick = [];
view([40 20]);