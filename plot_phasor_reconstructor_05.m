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


%% %%%%%%%%%%% Gráfica de reconstrucción animada
figure,
ax1=subplot(1,2,1);
hold on,
ax2=subplot(1,2,2);
plot(theta_i,rho_i,'b'),
hold on,
% xlim([min(theta_i),max(theta_i)])
% ylim([mini,maxi])
set(gcf ,'units','normalized')
set(gcf ,'position',[0.0702    0.5057    0.5958    0.3962])

rho_r=0;

nhar=5;
phasorcell={nhar-1,4};
for har=2:nhar
    rho_e=rhoiF(har).*exp(-1i*2*pi*f(har).*t);
    rho_real=real(rho_e);
    rho_imag=imag(rho_e);
    [angle,mag]=cart2pol(rho_real,rho_imag);
    phasorcell{har-1,1}=rho_real;
    phasorcell{har-1,2}=rho_imag;
    phasorcell{har-1,3}=angle;
    phasorcell{har-1,4}=mag;
end 
close all
figure,
hold on,
xlim([-50,50])
ylim([-50,50])
for i=1:length(rho_real)
    xp=0;yp=0;
    for har=2:nhar
    x=phasorcell{har-1,4}(i)*cos(phasorcell{har-1,3}(i));
    y=phasorcell{har-1,4}(i)*sin(phasorcell{har-1,3}(i));
    yc=yp+y; xc=xp+x;
    l(har-1)=plot([yp,yc],[xp,xc],'k');
    l2(har-1)=plot(yp+phasorcell{har-1,2},xp+phasorcell{har-1,1},'k');
    yp=yc; xp=xc;
    end
    drawnow
    delete(l)
    delete(l2)
end
return
figure,
hold on,



i




for har=2:128
    rho_e=rhoiF(har).*exp(-1i*2*pi*f(har).*t);
    rho_r=rho_r+rhoi_mag(har)*.2.*cos(2*pi*f(har).*t-rhoi_angle(har));
    rho_real=real(rho_e);
    rho_imag=imag(rho_e);
    [angle,mag]=cart2pol(rho_real,rho_imag);
    subplot(1,2,1)
    for i=1:length(rho_real)
    plot(rho_imag(i),rho_real(i),'.r');
    end
    for i=1:length(rho_real)
    x=mag(i)*cos(angle(i));
    y=mag(i)*sin(angle(i));
    l=quiver(0,0,y,x,'b');
    set(l,'AutoScaleFactor',1.01)
    drawnow
    delete(l)
    end
%     
%     for i=1:length(rho_real)
%     plot(rho_imag(i),rho_real(i),'.r');
%     subplot(1,2,1)
%     plot(rho_imag(i),rho_real(i),'.r');
%     [xa1,ya1] = ds2nfu(ax1,rho_imag(i),rho_real(i));
%     x=mag(i)*cos(angle(i));
%     y=mag(i)*sin(angle(i));
%     l=quiver(0,0,y,x,'b');
%     set(l,'AutoScaleFactor',1.01)
%     subplot(1,2,2)
%     plot(theta_i(i),rho_r(i),'.r')
%     [xa2,ya2] = ds2nfu(ax2,theta_i(i),rho_r(i));
%     xl = [xa1, xa2];
%     yl = [ya1,  ya2];
%     an=annotation('line',xl,yl);
%     drawnow
%     delete(an)
%     delete(l)
%     end
%     subplot(1,2,1)
%     xlim([min(theta_i),max(theta_i)])
%     ylim([mini,maxi])
%     plot(p_imag(i),p_real(i),'.r');
%     [xa1,ya1] = ds2nfu(ax1,p_imag(i),p_real(i));
%     x=mag(i)*cos(angle(i))/L;
%     y=mag(i)*sin(angle(i))/L;
%     l=quiver(0,0,y,x,'b');
%     set(l,'AutoScaleFactor',1.01)
%     subplot(1,2,2)
%     plot(theta_i(i),rho_r(i),'.r')
%     [xa2,ya2] = ds2nfu(ax2,t(i),recupera(i));
%     xl = [xa1, xa2];
%     yl = [ya1,  ya2];
%     an=annotation('line',xl,yl);
    drawnow
    
    
end



for har=2:128
    ax2=subplot(1,2,2);
    hold off
    plot(theta_i,rho_i,'b')
    hold on
    rho_r=rho_r+rhoi_mag(har).*2.*cos(2*pi*f(har).*t+rhoi_angle(har));
    if mod(har,40)==0
    for i=1:1:length(theta_i)
        % subplot(1,2,1)
        % plot(p_imag(i),p_real(i),'.r');
        % [xa1,ya1] = ds2nfu(ax1,p_imag(i),p_real(i));
        % x=mag(i)*cos(angle(i))/L;
        % y=mag(i)*sin(angle(i))/L;
        % l=quiver(0,0,y,x,'b');
        % set(l,'AutoScaleFactor',1.01)
        subplot(1,2,2)
        plot(theta_i(i),rho_r(i),'.r')
        % [xa2,ya2] = ds2nfu(ax2,t(i),recupera(i));
        % xl = [xa1, xa2];
        % yl = [ya1,  ya2];
        % an=annotation('line',xl,yl);
        drawnow
        % delete(an)
        % delete(l)
    end
    end
end




