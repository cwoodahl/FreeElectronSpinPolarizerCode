clear all; %close all;
%user inputs
beta = 0.1;
lambda =2.4e-6;
g1 = 0.475; phigm1 = 0; %interaction strenghts
g2 = 0.475; phigm2 = pi/2;
ns = 75; %pixel number

g1s = linspace(0,3,ns);
g2s = linspace(0,3,ns);
Ls = linspace(0,1,ns);

for c1 = 1:length(g1s)
    g1 = g1s(c1);
    g2 = g1; %commment out if not desired
    for c2 = 1:length(g2s)
        %g2 = g2s(c2);
        phipm = pi/2; Lqr = 1;
        %L = 0.5*Lqr;
        L = Ls(c2);


        if(phipm == 0)
            fac = 1;
        else
            fac = -1;
        end
        lambdag = beta*lambda;
        kz = 2*pi/lambdag;
        oz = 2*lambdag;
        z = linspace(-oz*2.5,oz*2.5,2000);
        sumT = 0; sumT2 = 0;
        ind = 16;
        for n = -ind:ind
            sum1 = besselj(n,2*g1)*exp(1i*n*(phigm1+phipm+kz*z)-1i*n^2*pi*L/Lqr);
            sumT = sumT + sum1;
        end
        sum2 = sin(2*2*g2*sin(kz*z+phigm2));
        %figure(1); hold on; plot(z,sum2,'LineWidth',2);

        dz = z(2)-z(1);
        norF = sum(sumT.*conj(sumT).*exp(-z.^2/2/oz^2));
        PsiD = sumT.*conj(sumT).*exp(-z.^2/2/oz^2)/norF;
        %figure(2); hold on; plot(z,PsiD,'LineWidth',2)

        sumT = fac*PsiD.*sum2;
        %figure(3); hold on; plot(z,sumT,'LineWidth',2)

        expect = sum(sumT);

        expects(c1,c2) = expect;
        %{
        figure; hold on;
        plot(z,sum2/max(sum2),'LineWidth',2);
        plot(z,PsiD/max(PsiD),'LineWidth',2);
        %}
    end
end

f = figure;
f.Position = [600 300 280 230]; 
imagesc(Ls,g1s,expects);
col = colorbar('eastoutside','XTickLabel',{'-0.7','0','0.7'},'XTick',[-0.7,0,0.7]);
colormap('jet');
ylabel(col,'\langle S_y \rangle','FontSize',12,'fontweight','normal');
%title('L_D = L_{QR}/2','FontSize',12,'fontweight','normal')
title('|g_{B1}| = |g_{B2}| \approx 0.5','FontSize',12,'fontweight','normal')
%ylabel('|g_{B1/2}|','FontSize',12);
%xlabel('L_D','FontSize',12);
ylabel('|g_{B1/2}|','FontSize',12);
xlabel('L_D/L_{QR}','FontSize',12);
ax = gca;
ax.FontSize = 12;
set(gca,'YDir','normal');
caxis([-0.7,0.7])

    
