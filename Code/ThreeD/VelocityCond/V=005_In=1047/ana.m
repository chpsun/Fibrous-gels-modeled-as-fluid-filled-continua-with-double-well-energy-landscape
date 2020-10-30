close all
clear all
clc

load('Results')
load('Parameters')

L = 10.47;
repout = 600;
%repin =  16500;

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



%TotalStep = repout * repin;
dL = -velocity * repout * time_Lstep ;



Time = 1:length(StressForward);

Strain_Compressive = Time /length(StressForward) * (dL) / L;

figure(1)
hold on
box on
set(gca,'FontSize',28);
plot(Strain_Compressive, -StressForward*kBT/v,'-b','linewidth',2)
plot(Strain_Compressive, -StressBackward(end:-1:1,:)*kBT/v,'-r','linewidth',2)
xlabel('Compressive Strain')
ylabel('Compressive Stress (Pa)')
legend('Compression','Stretching')


%%
[size1 size2] = size(StrainBackward);
x = linspace(0,1,size2);
t = linspace(0,10,size1);

figure(4)
hold on
box on
set(gca,'FontSize',28);
Matrix = [1 300 400 600];
for j = 1 :4
    i = Matrix(j);
    plot(x,1-StrainForward(i,:)/10.48,'linewidth',3)
end
xlabel('Reference Height X^*=X/H_0')
ylabel('Compressive Stretch \eta')
%legend('T=1s','T=201s','T=401s','T=601s','T=801s')
%T1 = title('$\dot{\gamma}=0.05 s^{-1}$');
%set(T1,'interpreter','latex');
axis([0 1 0 1])
legend('T=    0s','T=  90s','T=120s','T=180s')


%%

figure(6)
hold on
box on
set(gca,'FontSize',28);
Matrix = [1 300 401 600];
for j = 1 :4
    i = Matrix(j);
    plot(x,1-StrainBackward(i,:)/10.48,'linewidth',3)
end
xlabel('Reference Height X^*=X/H_0')
ylabel('Compressive Stretch \eta')
%legend('T=1s','T=201s','T=401s','T=601s','T=801s')
%T1 = title('$\dot{\gamma}=0.05 s^{-1}$');axis([0 1 0 1])
legend('T=    0s','T=  90s','T=120s','T=180s')