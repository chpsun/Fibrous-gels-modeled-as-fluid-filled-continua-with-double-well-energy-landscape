clear all
%%function pdex11
% kBT = 1.38*10^(-23) * 300;
% v = 1.7 * 10 ^ (-28);
% Ki = 0.1;
% G0 = 0.1 * 10 ^6;
% D = 5 * 10^(-9);
% coefa = 1.4345e+03;
% coefb = -12 * coefa;
% coefc = 36 * coefa;
% 
% etaG = v*G0/kBT;
% etac = v*coefc/kBT;
% etab = v* coefb/kBT;
% etaa = v* coefa/kBT;
rep = 10 ^ 3;
mul = 10 ^3;
m = 0;
x = linspace(0,1,10*10^1);
t = linspace(0,5,rep*mul);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
for i = 1 : rep
    u(i,:) = sol(mul*i,:,1);
end
% A surface plot is often a good way to study a solution.
% surf(x,t,u) 
% title('Numerical solution computed with 20 mesh points.')
% xlabel('Distance x')
% ylabel('Time t')

% A solution profile can also be illuminating.
% figure
% plot(x,u(end,:))
% title('Solution at t = 2')
% xlabel('Distance x')
% ylabel('u(x,2)')
% --------------------------------------------------------------
save('solution','u')
function [c,f,s] = pdex1pde(x,t,u,DuDx)
kBT = 1.38*10^(-23) * 300;
v = 1.7 * 10 ^ (-28);
Ki = 0.1;
G0 = 0.1 * 10 ^6;
D = 5 * 10^(-9);
% coefa = 794.4710;
% coefb = 8 * coefa;
% coefc = 7.25 * coefa;
% coefa =1.2*10^3;
% coefb =-12 * coefa;
% coefc = 36 * coefa;

coefa = 1.0 * 10^4;          % Material Property
coefb = -16 * coefa;
coefc = 64 * coefa;


c = 1/u^2 ;

f = (1-u)/u *DuDx * ...
    ( (1/(1-u)+1+2*Ki*u)-v*G0/kBT-v*coefc/kBT*u^(-2)...
    -2*v*coefb/kBT*(1/u-1)*u^(-2)-3*v*coefa/kBT*(1/u-1)^2*u^(-2));

%f = (1-u)/u *DuDx * ...
%    (1 * (1/(1-u)+1+2*Ki*u)-etaG-etac*u^(-2) - 2*etab*(1/u-1)*u^(-2)-3*etaa*(1/u-1)^2*u^(-2));

s = 0;
end
%% This is the initial condition
% --------------------------------------------------------------
function u0 = pdex1ic(x)


u0 = 0.7 + 2 /(1+exp(-40*x)) * (-0.3); 
end
%% This is the boundary condition
% --------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul-0.4;
ql = 0;
pr = 0;
qr = 1;
end