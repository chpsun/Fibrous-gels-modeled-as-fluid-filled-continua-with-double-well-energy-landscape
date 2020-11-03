%% For compression of gel
%% Goal is to plot constitutive laws

c2 = 6.4e5;     %% Units is Pa
c3 =-1.6e5;     %% Units is Pa
c4 = 1e4;       %% Units is Pa
G0 = 1e3;       %% I don't know the correct G0 

%% Now consider the pure one-dimensional compression 
phi = [0.05:0.001:1];    %% Solid volume fraction
lambda = (1/sqrt(3))*sqrt(phi.^(-2) + 2);    %% In pure one dimension 
%% lambda = 1 is the dry state, 10 is the fully swollen initial state. 
lambda1 = 1./phi;      %% This is the lambda_1 of PRSA paper

f4 = (c2/2)*(lambda-1).^2 + (c3/3)*(lambda-1).^3 + (c4/4)*(lambda-1).^4;
psi = f4 - G0*log(1./phi);
sigma = (c2*(lambda-1) + c3*(lambda-1).^2 + c4*(lambda-1).^3)./(3*lambda)./phi % need to multiply by 1/(3lambda phi)
%% Need to add the contribution of G0*log(detF) to the sigma above
psiold = (c2/2)*(lambda).^2 + (c3/3)*(lambda).^3 + (c4/4)*(lambda).^4;
strain = lambda1/max(lambda1) - 1;    %% Strain in expt

figure(1); 
plot(strain,psi,'k-','Linewidth',1.5);
hold on;
%%plot(strain,psiold,'r-','Linewidth',1.5);
xlabel('Expt. strain','Fontsize',16);

figure(2); 
plot(-strain,-sigma,'b-','Linewidth',1.5);
hold on;
xlabel('Expt. strain','Fontsize',16);
ylabel('Expt. stress (Pa)','Fontsize',16);

