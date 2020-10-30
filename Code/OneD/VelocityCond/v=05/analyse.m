% Here we assume a height of 100um.
% with D = 5 * 10^(-9)


clear all
close all

load('solutionCompress')
%%
[size1 size2] = size(u);
x = linspace(0,1,size2);
t = linspace(0,10,size1);

for i = 1 : size2-1
    eps(:,i) = (u(:,i+1)-u(:,i))/(x(i+1)-x(i));
end
eps(:,size2) = eps(:,size2-1);




%%
figure(1)
plot(x,eps(1:100:901,:))
xlabel('Reference Hight X^*=X/H_0')
ylabel('Stretch \lambda')

%%
figure(2)
plot(t,u(:,end))
xlabel('time t^*')
ylabel('Current Hight h^* = h/H_0')

%%
figure(3)
plot(x,u(1:100:901,:))
xlabel('Reference Hight X^*=X/H_0')
ylabel('Displacement u')

%%
figure(4)
hold on
box on
set(gca,'FontSize',28);
Matrix = [1 303 600 900];
for j = 1 :4
    i = Matrix(j);
    plot(x,1-eps(i,:)/10,'linewidth',3)
end
xlabel('Reference Height X^*=X/H_0')
ylabel('Compressive Stretch \eta')
%legend('T=1s','T=201s','T=401s','T=601s','T=801s')
%T1 = title('$\dot{\gamma}=0.05 s^{-1}$');
%set(T1,'interpreter','latex');

legend('T=      2s','T=  600s','T=1200s','T=1800s')


%% force extension
%  This is for compression
kBT = 1.38*10^(-23) * 300;
v = 1.7 * 10 ^ (-28);
Ki = 0.1;
G0 = 0.1 * 10 ^6;
D = 5 * 10^(-9);
% coefa = 794.4710;
% coefb = 8 * coefa;
% coefc = 7.25 * coefa;
coefa = 1.0 * 10^4;          % Material Property
coefb = -16 * coefa;
coefc = 64 * coefa;


DuDx = eps(:,end);

f_C = (log(1-1./DuDx) + 1./DuDx + Ki./DuDx.^2) ...
    -v*G0/kBT./DuDx + v*coefc/kBT* (DuDx-1)...
    +v*coefb/kBT*(DuDx-1).^2 + v*coefa/kBT*(DuDx-1).^3;

%% Stretching
u_C = u;
lambda_C = (10- u_C(:,end))/10;
load('solutionStrech')
[size1 size2] = size(u);
x = linspace(0,1,size2);
t = linspace(0,10,size1);
u_S = u;
lambda_S = (10-u_S(:,end))/10;
for i = 1 : size2-1
    eps_s(:,i) = (u(:,i+1)-u(:,i))/(x(i+1)-x(i));
end
eps_s(:,size2) = eps_s(:,size2-1);

DuDx = eps_s(:,end);

f_S = (log(1-1./DuDx) + 1./DuDx + Ki./DuDx.^2) ...
    -v*G0/kBT./DuDx + v*coefc/kBT* (DuDx-1)...
    +v*coefb/kBT*(DuDx-1).^2 + v*coefa/kBT*(DuDx-1).^3;

figure(5)
hold on
box on
set(gca,'FontSize',28);
SmfC = smooth(-f_C([1:25:end,end])*kBT/v);
SmfS = smooth(-f_S([1:25:end,end])*kBT/v);
plot(lambda_C([1:25:end,end]),SmfC,'linewidth',3);
plot(lambda_S([1:25:end,end]),SmfS,'linewidth',3);
xlabel('Compressive Strain')
ylabel('Compressive Stress (Pa)')
legend('Compression','Stretching')
%T1 = title('$\dot{\gamma}=0.005 s^{-1}$');
%set(T1,'interpreter','latex');
hold off

%%
 PB = 36;
xPB1 = x(1:PB);
u_CPB1 = u_C(end,1:PB);
xPB2 = x(PB:end);
u_CPB2 = u_C(end,PB:end);
%plot(xPB1,u_CPB1)
%plot(xPB2,u_CPB2)

figure(6)
hold on
box on
set(gca,'FontSize',28);
Matrix = [1 300 600 900];
for j = 1 :4
    i = Matrix(j);
    plot(x,1-eps_s(i,:)/10,'linewidth',3)
end
xlabel('Reference Height X^*=X/H_0')
ylabel('Compressive Stretch \eta')
%legend('T=1s','T=201s','T=401s','T=601s','T=801s')
%T1 = title('$\dot{\gamma}=0.05 s^{-1}$');
legend('T=      2s','T=  600s','T=1200s','T=1800s')
%set(T1,'interpreter','latex');


