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
load('solutionCompress');
rep = 10 ^ 3;
mul = 10 ^3;
m = 0;
x = linspace(0,1,10*10^1);
t = linspace(0,10,rep*mul);

ini = u(end,:);
clear u
%%
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
for i = 1 : rep
    u(i,:) = real(sol(mul*i,:,1));
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
save('solutionStrech','u')
save('solStrech','sol')
function [c,f,s] = pdex1pde(x,t,u,DuDx)
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


c = 1/(DuDx-1);

f = (log(1-1/DuDx) + 1/DuDx + Ki/DuDx^2) ...
    -v*G0/kBT/DuDx + v*coefc/kBT* (DuDx-1)...
    +v*coefb/kBT*(DuDx-1)^2 + v*coefa/kBT*(DuDx-1)^3;

%f = (1-u)/u *DuDx * ...
%    (1 * (1/(1-u)+1+2*Ki*u)-etaG-etac*u^(-2) - 2*etab*(1/u-1)*u^(-2)-3*etaa*(1/u-1)^2*u^(-2));

s = 0;
end
%% This is the initial condition
% --------------------------------------------------------------
function u0 = pdex1ic(x)


u0 = (x<=0.3535) * 9.443 * x + (x>0.3535) * (1.005*x^3 -3.547 * x^2 + 5.888 * x + 1.6555);
end
%% This is the boundary condition
% --------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur - (5 +0.5*t); % v^* = 500um/s
qr = 0;
end