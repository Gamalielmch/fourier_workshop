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
imshow(bw)

%% %%%%%%% Método Radius-Theta Centroide, perimetro 
stats=regionprops(bw,'Area','PixelIdxList','Centroid'); 
[~,ind]=max([stats.Area]);
%ind=1;
bw=bw*0;
bw(stats(ind).PixelIdxList)=1;
Cen=stats(ind).Centroid;
bws= bwboundaries(bw);
figure, hold on, axis ij
plot(Cen(1),Cen(2),'bs')
rho=zeros(1,length(bws{1})-1);
theta=zeros(1,length(bws{1})-1);
for i=1:length(bws{1})-1
plot(bws{1}(i,2),bws{1}(i,1),'r.')
h=plot([Cen(1) bws{1}(i,2)], [Cen(2) bws{1}(i,1)]);
difx=Cen(1)-bws{1}(i,2);
dify=Cen(2)-bws{1}(i,1);
[theta(i),rho(i)]=cart2pol(difx,dify);
drawnow
delete(h)
end

%% %%%%%%%Inter/extrapolation, perimeter
figure, subplot(1,2,1),polar(theta,rho)
np = 256; 
theta_i=linspace(-pi,pi,np);
rho_i=interp1(theta,rho,theta_i);
subplot(1,2,2), plot(theta,rho), hold on
plot(theta_i,rho_i,'r--')


