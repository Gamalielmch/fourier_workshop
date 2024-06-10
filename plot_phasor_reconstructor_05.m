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







%% %%%%%%%%%%% Gráfica de la señal sintetizada
% rho_r=0;
% rho_r=rho_r+rhoi_mag(1).*cos(2*pi*f(1).*t+rhoi_angle(1));
% subplot(1,2,1),title(['harmonics: ', num2str(0)])
% h=plot(theta_i,rho_r,'r--');
% subplot(1,2,2)
% h2=polar(theta_i,rho_r,'r');
% pause(2)
% delete(h)
% delete(h2)

% for i=2:n/2
% rho_r=rho_r+rhoi_mag(i)*2.*cos(2*pi*f(i).*t+rhoi_angle(i));
% subplot(1,2,1),title(['harmonics: ', num2str(i-1)])
% h=plot(theta_i,rho_r,'r--');
% subplot(1,2,2)
% h2=polar(theta_i,rho_r,'r');
% pause(2)
% delete(h)
% delete(h2)
% end

%% %%%%%%%%%%% Gráfica de reconstrucción animada
figure,
ax1=subplot(1,2,1);
hold on,
ax2=subplot(1,2,2);
plot(theta_i,rho_i,'b'),
hold on,
xlim([min(theta_i),max(theta_i)])
ylim([mini,maxi])
set(gcf ,'units','normalized')
set(gcf ,'position',[0.0702    0.5057    0.5958    0.3962])

rho_r=0;
rho_r=rho_r+rhoi_mag(1).*cos(2*pi*f(1).*t+rhoi_angle(1));

harmonic=2;
p=rhoiF(harmonic).*exp(-1i*2*pi*f(harmonic).*t);
p_real=real(p);
p_imag=imag(p);
[angle,mag]=cart2pol(p_real,p_imag);
rho_r=rho_r+rhoi_mag(2)*cos(2*pi*f(2).*t-rhoi_angle(2));
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




