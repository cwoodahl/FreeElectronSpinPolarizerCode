%clear all; %close all;
beta = 0.1; %electron relativistic fraction
lambda =2.4e-6; %laser wavelength

g1 = 0.475; phigm1 = 0; %first interaction stage strength and phase
g2 = 0.475; phigm2 = pi/2; %second interaction stage strength and phase
phipm = pi; %spin state phase
Lqr = 1; %talbot revival length
L = 0.5*Lqr; %L given in terms of talbot revival length

if(phipm == 0)
    fac = 1;
    colSpin = [0,0.5,0.8];
    lsty = '-';
    lW = 1.8;
else
    fac = -1;
    colSpin = [0.7,0,0.7];
    lsty = '-.';
    lW = 1.6;
end
lambdag = beta*lambda;
kz = 2*pi/lambdag;
oz = 2*lambdag;
z = linspace(-oz*2.5,oz*2.5,10000);
sumT = 0; sumT2 = 0;
ind = 1200;
for n = -ind:ind
    sum1 = 1i^(n)*besselj(n,2*g1)*exp(1i*n*(phigm1+phipm+kz*z)-1i*n^2*pi*L/Lqr);
    sumT = sumT + sum1;
end
sum2 = sin(2*2*g2*sin(kz*z+phigm2));
%f0=figure(20); hold on; plot(z,sum2,'LineWidth',1.6);
%f0.Position = [600 300 235 250]; 

dz = z(2)-z(1);
norF = sum(sumT.*conj(sumT).*exp(-z.^2/2/oz^2));
PsiD = sumT.*conj(sumT).*exp(-z.^2/2/oz^2)/norF;
f1 = figure(22); hold on; plot(z/(1e-6),PsiD,'LineStyle',lsty,'LineWidth',lW,'color',colSpin)
f1.Position = [600 300 310 164]; 
xlabel('z^, (\mum)','fontsize',12.5)
ylabel('Intensity (arb. units)')
title('Probability Density','fontsize',12.5)
ax = gca;
ax.FontSize = 12;

sumT = fac*PsiD.*sum2;
f2=figure(23); hold on; plot(z/(1e-6),sumT,'LineWidth',lW,'color',colSpin,'LineStyle',lsty)
f2.Position = [600 300 310 170]; 
xlabel('z^, (\mum)','fontsize',12.5)
ylabel('Intensity (arb. units)')
title('\langle S_y \rangle Density, \phi_g = \pi','fontsize',12.5)
ax = gca;
ax.FontSize = 12;


expect = sum(sumT)

%figure; hold on;
%plot(z,sum2/max(sum2),'LineWidth',1.6);
%plot(z,PsiD/max(PsiD),'LineWidth',1.6);
            
