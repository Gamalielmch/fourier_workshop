clear all
close all
I=imread('goku.png');
I=imresize(I,0.3);
bw=imbinarize(I);
bw = imfill(bw,'holes');
bw = bwmorph(bw,'majority');
bw = bwmorph(bw,'majority');
bw = bwmorph(bw,'majority');
bw=bwperim(bw);


bw=padarray(bw,[2,2],0);
size_I=size(I);
imshow(bw)
bw=bwperim(bw);

imshow(bw)
[yo,xo]=find(bw==1,1);
lon=sum(double(bw(:)));
cadena=zeros(1,lon);
%bwr=bw*0;
% for i=1:lon
%     bwr(bws(i,1),bws(i,2))=1;
%     imshow(bwr)
%     drawnow
% end
xi=xo;
yi=yo;
le=sqrt(2);
for i=1:lon-1
    try
    pos=[yi,xi+1;yi-1,xi+1;yi-1,xi;yi-1,xi-1;yi,xi-1;yi+1,xi-1;yi+1,xi;yi+1,xi+1];
    neigh=[bw(pos(1,1),pos(1,2)),bw(pos(2,1),pos(2,2)),bw(pos(3,1),pos(3,2)),...
        bw(pos(4,1),pos(4,2)),bw(pos(5,1),pos(5,2)),bw(pos(6,1),pos(6,2)),...
        bw(pos(7,1),pos(7,2)),bw(pos(8,1),pos(8,2))];
    [~,mino]=max(neigh);
    cadena(i)= mino(end)-1;
    bw(yi,xi)=0;
    xi=pos(mino,2);
    yi=pos(mino,1);
    imshow(bw), drawnow
    catch
        i
    end
end



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
