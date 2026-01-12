%% constants
clear all;
c=3e8; h=6.626e-34; hbar=h/2/pi; me=9.11e-31; q=1.6e-19; mub = 9.274e-24;

col = 'red';
%% user inputs
fac=3.85;
L = 12.7e-6*fac; %interaction propagation length
Ey0 = 385e6/fac; %field intensity 
lambda = 2.4e-6; %light wavelength
beta =0.1; %electron relativsitic fraction
sigmaz = 250e-9; %rms electron width
tSteps = 150; %number of time steps to split into
phi0 = pi; %interaction phase
Lr = 0.5; %portion of Lqr drifted before second interaction
g1 = 0.475; %first interaction state coupling factor
czp = 0; %z spin up state coeffcient
czm = 1; %z spin down state coeffcient
L2 = 0; %portion of Lqr between end of interaction region and measurement
phi1 = 0; %first stage interaction phase

%% calculations from inputs
omega = 2*pi*c/lambda;
Ay0 = Ey0/omega;
v0 = beta*c; vo = v0;
tp = L/v0;
kz = 2*pi/lambda/beta;
dtp = tp/tSteps;
gamma = 1/sqrt(1-beta^2);
Lqr = beta^3*gamma^3*me*c*lambda^2/h;
L2 = L2*Lqr;
dtd2 = L2/v0;



%% initialize electron wavefunction
Nz = 3048;
zps = linspace(-3.4*sigmaz,3.4*sigmaz,Nz);
dzps = zps(2)-zps(1); dz = dzps;

sumV1 = 0; sumV2 = 0;
for n = -3:3
    sumV1 = sumV1 + besselj(n,2*g1).*exp(1i*n*kz*zps+1i*n*phi1)...
        .*exp(-1i*n^2*pi*Lr);
    sumV2 = sumV2 + besselj(n,2*g1).*exp(1i*n*kz*zps+1i*n*(phi1+pi))...
        .*exp(-1i*n^2*pi*Lr);
end

PsiU = exp(-zps.^2/4/sigmaz^2).*sumV1*czp;
PsiD = exp(-zps.^2/4/sigmaz^2).*sumV2*czm;
%figure;plot(zps,abs(PsiU).^2-abs(PsiD).^2,'LineWidth',1.5);

ks = linspace(-24*kz,24*kz,1000);
PsiU_k = zeros(1,length(ks)); PsiD_k = zeros(1,length(ks));
for ctr = 1:length(ks);
    PsiU_k(ctr) = sum(PsiU.*exp(1i*ks(ctr).*zps));
end

for ctr = 1:length(ks);
    PsiD_k(ctr) = sum(PsiD.*exp(1i*ks(ctr).*zps));
end

%figure; plot(ks,abs(Psi_k).^2,'LineWidth',1.5)

%% Precalc
    %set up operators in space 
    Ay = Ay0*exp(1i*kz*zps+1i*phi0)+Ay0*exp(-1i*kz*zps-1i*phi0);
    Bx = Ay0*kz*-1i*exp(1i*kz*zps+1i*phi0)+Ay0*1i*kz*exp(-1i*kz*zps-1i*phi0);

    g = mub/hbar*dtp/2*Bx;
    %Ay = 0;
    UA = exp(-1i*dtp/2/hbar*q^2*Ay.*conj(Ay)/2/me);
    
    phase = exp(-1i * ks.' * zps);
    phase2 = exp(1i * (ks.' * zps));
    Eq = hbar^2.*ks.^2/2/gamma^3/me;
    UD = exp(-1i*Eq*dtp/hbar);
    UD0 = UD;
    %UD = 1;

%% calculate exp(-1i*N*dt/hbar)

for ctr2 = 1:tSteps
    
    PsiUo = PsiU; PsiDo = PsiD;
    PsiU = (1-g.^2/2).*PsiUo-(1i*g).*PsiDo;
    PsiD = (1-g.^2/2).*PsiDo-(1i*g).*PsiUo;
    
    PsiU = UA.*PsiU;
    PsiD = UA.*PsiD;
    
    %now go to k space
    %for ctr = 1:length(ks)
    %    PsiU_k(ctr) = sum(PsiU.*exp(1i*ks(ctr).*zps));
    %    PsiD_k(ctr) = sum(PsiD.*exp(1i*ks(ctr).*zps));
    %end
    
    %
    PsiU_k2 = phase2 * PsiU.';  % (1000×2048)*(2048×1) = (1000×1)
    PsiU_k = PsiU_k2.';
    PsiD_k2 = phase2 * PsiD.';  % (1000×2048)*(2048×1) = (1000×1)
    PsiD_k = PsiD_k2.';
    
    
    %quadratic phase accumulation factor in k-space
    PsiU_k = UD.*PsiU_k;
    PsiD_k = UD.*PsiD_k;
    
    %now go back to space space
    %for ctr = 1:length(zps)
    %    PsiU(ctr) = sum(PsiU_k.*exp(-1i.*ks*zps(ctr)));
    %    PsiD(ctr) = sum(PsiD_k.*exp(-1i.*ks*zps(ctr)));
    %end
    
    % Build phase matrix and perform vectorized inverse transform
    PsiU = PsiU_k * phase;  % row (1×Nk) × (Nk×Nz) = (1×Nz)
    PsiD = PsiD_k * phase;


    
    PsiUo = PsiU; PsiDo = PsiD;
    PsiU = (1-g.^2/2).*PsiUo-(1i*g).*PsiDo;
    PsiD = (1-g.^2/2).*PsiDo-(1i*g).*PsiUo;
    
    PsiU = UA.*PsiU;
    PsiD = UA.*PsiD;
    
    norm0 = sqrt(sum(abs(PsiU).^2+abs(PsiD).^2));
    PsiU = PsiU/norm0;
    PsiD = PsiD/norm0;
end

%% free space drift after interaction region
PsiU_k2 = phase2 * PsiU.';  % (1000×2048)*(2048×1) = (1000×1)
PsiU_k = PsiU_k2.';
PsiD_k2 = phase2 * PsiD.';  % (1000×2048)*(2048×1) = (1000×1)
PsiD_k = PsiD_k2.';

UD = exp(-1i*Eq*dtd2/hbar);
PsiU_k = UD.*PsiU_k;
PsiD_k = UD.*PsiD_k;

PsiU = PsiU_k * phase;  % row (1×Nk) × (Nk×Nz) = (1×Nz)
PsiD = PsiD_k * phase;

norm0 = sqrt(sum(abs(PsiU).^2+abs(PsiD).^2));
PsiU = PsiU/norm0;
PsiD = PsiD/norm0;

%% plot functions
figure(11); subplot(3,1,1); hold on;
plot(zps*1e9,(abs(PsiU).^2+abs(PsiD).^2),'LineWidth',2,'color',col);
title('Spatial Probability Density')
xlabel('zp (nm)')

subplot(3,1,2); hold on;
plot(ks/kz,(abs(PsiU_k).^2+abs(PsiD_k).^2)/norm0^2/max((abs(PsiU_k).^2+abs(PsiD_k).^2)/norm0^2),'LineWidth',2,'color','black');
title('Momentum Offset Density')
xlabel('k')

subplot(3,1,3); hold on;
Sy = -1i*PsiD.*conj(PsiU)+1i*conj(PsiD).*PsiU;
expect = sum(Sy)
plot(zps*1e9,Sy/max(Sy),'LineWidth',2,'color',col);
title(['Y-Spin Expectation Density, P_y= ',num2str(round(expect*100,2)),'%'])
xlabel('zp (nm)')


    



