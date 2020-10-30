% Here we assume a height of 100um.
% with D = 5 * 10^(-9)

clear all
close all
load('solution')
%%
kBT = 1.38*10^(-23) * 300;
nu = 1.7 * 10 ^ (-28);
Ki = 0.1;
G0 = 0.1 * 10 ^6;
D = 5 * 10^(-9);
% coefa = 794.4710;
% coefb = 8 * coefa;
% coefc = 7.25 * coefa;
coefa = 1.0 * 10^4;          % Material Property
coefb = -16 * coefa;
coefc = 64 * coefa;
[size1 size2] = size(u);
x = linspace(0,1,size2);
t = linspace(0,1000,size1);

for i = 1 : size1
    H(i) = sum(1./u(i,:))/size2;
end

for i = 1 : size1-1
    v(i) = (H(i+1)-H(i))/(t(i+1)-t(i));
end
v(size1) = v(size1-1);
%%
figure(1)
plot(t,H)
xlabel('time t^*')
ylabel('Current Height h^* = h/H_0')
%%
figure(2)
plot(t,v/10,'linewidth',3)
box on
set(gca,'FontSize',28);
xlabel('time (s)')
ylabel('Strain rate (s^{-1})')
%title('S = 2.3\times 10^6 Pa');
%%
figure(3)
hold on
box on
set(gca,'FontSize',28);
plot(x(1,end:-1:1),1-1./u(1:200:801,:)/10,'linewidth',3)
xlabel('Reference Height X^*=X/H_0')
ylabel('Stretch \lambda')
%legend('T=1s','T=201s','T=401s','T=601s','T=801s')
%title('S = 2.3\times 10^6 Pa');

DuDx = 1/0.4; %1/0.4
((log(1-1./DuDx) + 1./DuDx + Ki./DuDx.^2) ...
    -nu*G0/kBT./DuDx + nu*coefc/kBT* (DuDx-1)...
    +nu*coefb/kBT*(DuDx-1).^2 + nu*coefa/kBT*(DuDx-1).^3)*kBT/nu

%% 4
figure(4)
plot(t * 2,v/20,'linewidth',3) % reference height H = 10; T*=2 s;
box on
set(gca,'FontSize',28);
xlabel('time (s)')
ylabel('Strain rate (s^{-1})')
xticks([0,1000,2000])
yticks([-0.010,-0.006,-0.002,0])
axis([0 2000 -8*10^(-3) 0]);
%
axes('Position',[0.37 0.3 0.5 0.5])
box on

plot(x(1,end:-1:1),1-1./u([1,201,401,601,801],:)/10,'linewidth',3)
xticks([0,0.2,0.4,0.6,0.8,1.0])
yticks([0,0.2,0.4,0.6,0.8])

xlabel('Reference height X^*=X/H_0','FontSize',22)
ylabel('Compressive stretch \lambda','FontSize',22)
set(gca,'FontSize',20);

Str = ['T =   2s'; 'T = 400s'; 'T = 800s'; 'T =1200s'; 'T =1600s'];
annotation('textbox',[.4 .3 .08 .05],'String',Str(1,:),'FitBoxToText','on','LineStyle','none','FontSize',14)
annotation('textbox',[.4 .46 .08 .05],'String',Str(2,:),'FitBoxToText','on','LineStyle','none','FontSize',14)
annotation('textbox',[.4 .65 .08 .05],'String',Str(3,:),'FitBoxToText','on','LineStyle','none','FontSize',14)
annotation('textbox',[.4 .71 .08 .05],'String',Str(4,:),'FitBoxToText','on','LineStyle','none','FontSize',14)
annotation('textbox',[.4 .752 .08 .05],'String',Str(5,:),'FitBoxToText','on','LineStyle','none','FontSize',14)







