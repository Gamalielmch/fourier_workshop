%% %%%%%%% Limpiar
clear all,
close all,
clc
%% %%%%%%% cargar imagen y binarizar imagen
I=imread('mdck.png');
bw=imbinarize(I(:,:,2),0.4);
bw = bwmorph(bw,'majority');
bw = bwmorph(bw,'majority');
bw = imfill(bw,'holes');


%% %%%%%%% Método Radius-Theta Centroide, perimeter
stats=regionprops(bw,'Area','PixelIdxList','Centroid');
[~,ind]=max([stats.Area]);
%ind=1;
bw=bw*0;
bw(stats(ind).PixelIdxList)=1;
Cen=stats(ind).Centroid;
bws= bwboundaries(bw);
rho=zeros(1,length(bws{1})-1);
theta=zeros(1,length(bws{1})-1);
for i=1:length(bws{1})-1
    difx=Cen(1)-bws{1}(i,2);
    dify=Cen(2)-bws{1}(i,1);
    [theta(i),rho(i)]=cart2pol(difx,dify);
end


%% %%%%%%%Inter/extrapolation, perimeter
np = 256;
theta_i=linspace(-pi,pi,np);
rho_i=interp1(theta,rho,theta_i,'linear','extrap');
rho_i(:)=rho_i(:)-min(rho_i(:));
rho_i=rho_i/max(rho_i(:));
rho_i=rho_i-mean(rho_i);
mini=min(rho_i(:));
maxi=max(rho_i(:));


%% %%%%%%%%%%% Parametros para la transformada
fs=length(theta_i);
T = 1/fs;             % Periodo de muestreo
L=length(theta_i);         % Longitud de la señal
t = (0:L-1)*T;         % Vector del tiempo
n = 2^nextpow2(L);     % Puntos de la tranformada
f = fs*(0:(n/2))/n;    % Vector de frecuencias

%% %%%%%%%%%%% Transformada de Fourier
rhoiF = fft(rho_i,n);        % Transformada
rhoi_mag = abs(rhoiF/L);
rhoi_angle = angle(rhoiF/L);


%% %%%%%%%%%%% Cálculo de fasores
nhar=5;
phasorcell={nhar-1,4};
rho_r=0;
rho_r=rho_r+rhoi_mag(1).*cos(2*pi*f(1).*t+rhoi_angle(1));

for har=2:nhar
    rho_e=rhoiF(har).*exp(1i*2*pi*f(har).*t)./128;
    rho_real=real(rho_e);
    rho_imag=imag(rho_e);
    rho_r=rho_r+rhoi_mag(har).*2.*cos(2*pi*f(har).*t+rhoi_angle(har));
    [angle,mag]=cart2pol(rho_real,rho_imag);
    phasorcell{har-1,1}=rho_real;
    phasorcell{har-1,2}=rho_imag;
    phasorcell{har-1,3}=angle;
    phasorcell{har-1,4}=mag;
end

%% %%%%%%%%%%% Gráfica de fasores

figure('Color',[1 1 1])
set(gcf ,'units','normalized')
set(gcf ,'position',[0.0339    0.3952    0.9155    0.5067])
ax1=subplot(1,2,1); hold on,
xlim([-1,1]), ylim([-0.3,0.8])
ax2=subplot(1,2,2);
hold on, ylim([-0.3,0.8])
plot(theta_i,rho_i)

for i=1:length(rho_real)
    xp=0;yp=0;
    subplot(1,2,1)
    for har=2:nhar
        x=phasorcell{har-1,4}(i)*cos(phasorcell{har-1,3}(i));
        y=phasorcell{har-1,4}(i)*sin(phasorcell{har-1,3}(i));
        yc=yp+y; xc=xp+x;
        
        l(har-1)=plot([yp,yc],[xp,xc],'k');
        l2(har-1)=plot(yp+phasorcell{har-1,2},xp+phasorcell{har-1,1},'k');
        yp=yc; xp=xc;
    end
    subplot(1,2,2)
    plot(theta_i(i),rho_r(i),'r.')
    [xa1,ya1] = ds2nfu(ax1, yc,xc);
    [xa2,ya2] = ds2nfu(ax2,theta_i(i),rho_r(i));
    xl = [xa1, xa2];
    yl = [ya1,  ya2];
    an=annotation('line',xl,yl);
    drawnow
    pause(0.1)
    delete(l)
    delete(l2)
    delete( an)
end


