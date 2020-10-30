% Several Targets in this program
% Start from a uniform configuration; Once reach a critical, a phase
% boundry will appear
% Backward Process start immediately after the forward ends

clear all
clc

T1 = clock;
N = 201;                     %Number of segments
L = 1;                      %Length of the sample
dx = L/(N-1);
dt = dx * 0.001;

%% Parameters

velocity1 =0.1 % 0.001;
R = L * 10;                 % Radius of the sample

kBT = 1.38*10^(-23) * 300;   % Temperature
v = 1.7 * 10 ^ (-28);        % Volume of liquid per atom
Ki = 0.1;
G0 = 0.1 * 10 ^6;            % A shear constant
D = 5 * 10^(-9);
coefa = 1.0 * 10^4;          % Material Property
coefb = -16 * coefa;
coefc = 64 * coefa;
coeff = 6;


vis = 0.0;

etaG = v*G0/kBT;
etac = v*coefc/kBT;
etab = v* coefb/kBT;
etaa = v* coefa/kBT;

% T_R = a (Du-1)^3 + b(Du-1)^2 + c(Du-1)

%% Matrix A
% A = zeros(N-4,N-4);
% for i = 1 : N-5
%     A(i,i+1) = 1;
%     A(i+1,i) = -1;
% end
% 
% A(1,1) = -1/2;
% A(N-4,N-4) = +1/2;
% 
% A_inv = inv(A);

%% initialization

X = linspace(0,L,N);
u = zeros(1,N);
u_new = u;
u_dot = zeros(1,N);

Du = u;
DDu = u;
%DDDu = u;

phi = u;
Dphi = u;
%DDphi = u;

mu = u;
Dmu = u;
%DDmu = u;
m = u;
%Dm = u;

S = 0; % S is dimensionless. vS/kBT

%% IC 
% a good combination is 9.53 for high strain phase and 3 for low 
% for v = 0.1 use 10.3
% u = X * 10;
N_PB = 195;
u(1:end) = X(1:end) * 10.40;
%  u(N_PB:end) = X(N_PB:end) * 3;
%  u(N_PB:end) = u(N_PB:end) + u(N_PB-1)-X(N_PB-1)*3;


repout = 300 ;%17 and imaginary
repin =  25000; %repin * velocity1 = 2500
repout_eq = 0; % 300;
time_Lstep = repin * dt;

Strain = zeros(repout,N);
Stress = zeros(repout,1);

%% Iteration Part

%% Forward
velocity = 0;
velocity = velocity1;
for time = 1 : repout
    for instep = 1 : repin
%% find the derivatives

for i = 2 : N-1
    
    Du(i) = ( u(i+1)-u(i) )/(dx);
    DDu(i) = ( u(i+1) - 2 * u(i) + u(i-1) ) / dx^2;
    Du_dot(i) = ( u_dot(i)-  u_dot(i) )/( dx);
%    DDDu(i) = ( u(i+2) - 2 * u(i+1) + 2 * u(i-1) - u(i-2) ) / dx^3;
end
    Du(1) = (u(2) - u(1))/dx;
%    Du(2) = ( u(2+1)-u(2-1) )/(2 * dx);
%    Du(N-1) = ( u(N)-u(N-2) )/(2 * dx);
    Du(N) = (u(N) - u(N-1))/dx;
    
%% Update S
Du_end = (u_new(N) - u_new(N-1))/dx;
phi_end = Du_end^(-1);

for i = 1 : N
    phi(i) = 1/Du(i);
    mu_S(i) = ( log(1-phi(i)) + phi(i) + Ki * phi(i) ^2) ...
             -v*G0/kBT/Du(i) + v*coefc/kBT* (Du(i)-1)...
             +v*coefb/kBT*(Du(i)-1)^2 + v*coefa/kBT*(Du(i)-1)^3 + vis * u_dot(i);
    m(i) = (1-phi(i))/phi(i);
end
         
         
S = (-velocity * pi * R^2/2  + coeff * sum( m.*mu_S) )/ sum(m);

    
%  phi, mu    
for i = 2 : N-1
%    phi(i) = 1/Du(i);
    Dphi(i) = -1/Du(i)^2 * DDu(i);
%    DDphi(i) = -1/Du(i)^2 * DDDu(i) + 2 /Du(i)^3 * DDu(i)^2;
    
    mu(i) = mu_S(i) - S;
    Dmu(i) = ((-1/(1-phi(i)) + 1 + 2 * Ki * phi(i) ) - v*G0/kBT ) * Dphi(i) + ...
             (v*coefc/kBT + 2 * v*coefb/kBT*(Du(i)-1) + 3 * v*coefa/kBT*(Du(i)-1)^2 ) * DDu(i) + vis * Du_dot(i);
    
%    DDmu(i) = ((-1/(1-phi(i)) + 1 + 2 * Ki * phi(i) ) - v*G0/kBT ) * DDphi(i) + ...
%             (1/(1-phi(i))^2 + 2 * Ki ) *Dphi(i)^2 +...
%             (v*coefc/kBT + 2 * v*coefb/kBT*(Du(i)-1) + 3 * v*coefa/kBT*(Du(i)-1)^2 ) * DDDu(i) + ...
%             ( 2 * v*coefb/kBT + 3 * v*coefa/kBT*(Du(i)-1) ) * DDu(i)^2;
             
%    m(i) = (1-phi(i))/phi(i);
%    Dm(i) = -phi(i)^(-2) * Dphi(i);

end

%% Update u
F = coeff * m(1) * mu(1);
for i = 2 : N-1
    F = F + coeff * m(i) * mu(i);
    u_new(i) = u(i) + dt * ( m(i) * Dmu(i) - 2/(pi*R^2) * dx * F    );
    u_dot(i) =  ( m(i) * Dmu(i) - 2/(pi*R^2) * dx * F    );
end

u_new(1) = 0;
u_new(N) = u(N) - velocity * dt;
u_dot(N) = -velocity;



u = u_new;
    end
Strain(time,:) = Du;
Stress(time,:) = S;
end

StrainForward = Strain;
StressForward = Stress;
clear Strain
clear Stress

%% Backward
velocity = - velocity;

for time = 1 : repout
    velocity = -velocity1 * (time>10) + (velocity1- 2*velocity1/10 * time) * (time<11);
    for instep = 1 : repin
%% find the derivatives

for i = 2 : N-1
    
    Du(i) = ( u(i+1)-u(i) )/(dx);
    DDu(i) = ( u(i+1) - 2 * u(i) + u(i-1) ) / dx^2;
    Du_dot(i) = ( u_dot(i)-  u_dot(i) )/( dx);
%    DDDu(i) = ( u(i+2) - 2 * u(i+1) + 2 * u(i-1) - u(i-2) ) / dx^3;
end
    Du(1) = (u(2) - u(1))/dx;
%    Du(2) = ( u(2+1)-u(2-1) )/(2 * dx);
%    Du(N-1) = ( u(N)-u(N-2) )/(2 * dx);
    Du(N) = (u(N) - u(N-1))/dx;
    
%% Update S
Du_end = (u_new(N) - u_new(N-1))/dx;
phi_end = Du_end^(-1);

for i = 1 : N
    phi(i) = 1/Du(i);
    mu_S(i) = ( log(1-phi(i)) + phi(i) + Ki * phi(i) ^2) ...
             -v*G0/kBT/Du(i) + v*coefc/kBT* (Du(i)-1)...
             +v*coefb/kBT*(Du(i)-1)^2 + v*coefa/kBT*(Du(i)-1)^3 + vis * u_dot(i);
    m(i) = (1-phi(i))/phi(i);
end
         
         
S = (-velocity * pi * R^2/2  + coeff * sum( m.*mu_S) )/ sum(m);

    
%  phi, mu    
for i = 2 : N-1
%    phi(i) = 1/Du(i);
    Dphi(i) = -1/Du(i)^2 * DDu(i);
%    DDphi(i) = -1/Du(i)^2 * DDDu(i) + 2 /Du(i)^3 * DDu(i)^2;
    
    mu(i) = mu_S(i) - S;
    Dmu(i) = ((-1/(1-phi(i)) + 1 + 2 * Ki * phi(i) ) - v*G0/kBT ) * Dphi(i) + ...
             (v*coefc/kBT + 2 * v*coefb/kBT*(Du(i)-1) + 3 * v*coefa/kBT*(Du(i)-1)^2 ) * DDu(i) + vis * Du_dot(i);
    
%    DDmu(i) = ((-1/(1-phi(i)) + 1 + 2 * Ki * phi(i) ) - v*G0/kBT ) * DDphi(i) + ...
%             (1/(1-phi(i))^2 + 2 * Ki ) *Dphi(i)^2 +...
%             (v*coefc/kBT + 2 * v*coefb/kBT*(Du(i)-1) + 3 * v*coefa/kBT*(Du(i)-1)^2 ) * DDDu(i) + ...
%             ( 2 * v*coefb/kBT + 3 * v*coefa/kBT*(Du(i)-1) ) * DDu(i)^2;
             
%    m(i) = (1-phi(i))/phi(i);
%    Dm(i) = -phi(i)^(-2) * Dphi(i);

end

%% Update u
F = coeff * m(1) * mu(1);
for i = 2 : N-1
    F = F + coeff * m(i) * mu(i);
    u_new(i) = u(i) + dt * ( m(i) * Dmu(i) - 2/(pi*R^2) * dx * F    );
    u_dot(i) =  ( m(i) * Dmu(i) - 2/(pi*R^2) * dx * F    );
end

u_new(1) = 0;
u_new(N) = u(N) - velocity * dt;
u_dot(N) = -velocity;



u = u_new;
    end
Strain(time,:) = Du;
Stress(time,:) = S;
end

StrainBackward = Strain;
StressBackward = Stress;

%%
save('Results','StrainForward','StressForward','StrainBackward','StressBackward');
save('Parameters', 'velocity','coefa','coefb','coefc','time_Lstep','dx')

T2 = clock;
etime(T2,T1)