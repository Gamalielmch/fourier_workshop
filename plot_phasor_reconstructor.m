clear all
close all
I=imread('mdck.png');
bw=imbinarize(I(:,:,2),0.4);
bw = bwmorph(bw,'majority');
bw = bwmorph(bw,'majority');
bw = imfill(bw,'holes');
%imshow(bw)
stats=regionprops(bw,'Area','PixelIdxList','Centroid'); 
[~,ind]=max([stats.Area]);
%ind=1;
bw=bw*0;
bw(stats(ind).PixelIdxList)=1;
Cen=stats(ind).Centroid;
bws= bwboundaries(bw);
%figure, hold on, axis ij
%plot(Cen(1),Cen(2),'bs')
rho=zeros(1,length(bws{1}));
theta=zeros(1,length(bws{1}));
for i=1:length(bws{1})
%plot(bws{1}(i,2),bws{1}(i,1),'r.')
%h=plot([Cen(1) bws{1}(i,2)], [Cen(2) bws{1}(i,1)]);
difx=Cen(1)-bws{1}(i,2);
dify=Cen(2)-bws{1}(i,1);
rho(i)=sqrt(difx^2+dify^2);
theta(i)=2*(atan(dify/difx));
%drawnow
%delete(h)
end
theta=unwrap(theta,0.9)/2;
%figure, polar(theta,rho)
np = 256; 
theta_i=linspace(0,2*pi,np);
rho_i=interp1(theta,rho,theta_i,'linear','extrap');
factor=max(rho_i);
rho_i=rho_i/factor;
%%%%%%%%%%%%% Parametros para la transformada 
fs=length(theta_i);
T = 1/fs;             % Periodo de muestreo   
L=length(theta_i);         % Longitud de la señal
t = (0:L-1)*T;         % Vector del tiempo
n = 2^nextpow2(L);     % Puntos de la tranformada 
f = fs*(0:(n/2))/n;    % Vector de frecuencias
%%%%%%%%%%%%% Transformada de Fourier señal uno
rhoiF = fft(rho_i,n);        % Transformada
rhoi_mag = abs(rhoiF/L); 
rhoi_angle = angle(rhoiF/L); 
% mag=sqrt(p_real.^2+p_imag.^2);
% angle=(atan(p_imag./p_real))*2;
% angle=unwrap(angle,0.9)/2;


figure,subplot(1,2,1),
plot(theta_i,rho_i,'b'), 
hold on
subplot(1,2,2)
polar(theta_i,rho_i,'b');
hold on

rho_r=0;
rho_r=rho_r+rhoi_mag(1).*cos(2*pi*f(1).*t+rhoi_angle(1));
h=plot(theta_i,rho_r,'k--');
title(num2str(i))
pause(0.5)
delete(h)

for i=2:n/2
rho_r=rho_r+rhoi_mag(i)*2.*cos(2*pi*f(i).*t+rhoi_angle(i));
subplot(1,2,1)
h=plot(theta_i,rho_r,'k--');
hold on
subplot(1,2,2)
h2=polar(theta_i,rho_r,'k');
title(num2str(i))
pause(2)
delete(h)
delete(h2)
end
