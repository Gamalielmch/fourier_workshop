clear all
close all
I=imread('mdck.png');
bw=imbinarize(I(:,:,2),0.4);
bw = bwmorph(bw,'majority');
bw = bwmorph(bw,'majority');
bw = imfill(bw,'holes');
imshow(bw)
stats=regionprops(bw,'Area','PixelIdxList','Centroid'); 
[~,ind]=max([stats.Area]);
%ind=1;
bw=bw*0;
bw(stats(ind).PixelIdxList)=1;
Cen=stats(ind).Centroid;
bws= bwboundaries(bw);
figure, hold on, axis ij
plot(Cen(1),Cen(2),'bs')
rho=zeros(1,length(bws{1}));
theta=zeros(1,length(bws{1}));
for i=1:length(bws{1})
plot(bws{1}(i,2),bws{1}(i,1),'r.')
h=plot([Cen(1) bws{1}(i,2)], [Cen(2) bws{1}(i,1)]);
difx=Cen(1)-bws{1}(i,2);
dify=Cen(2)-bws{1}(i,1);
rho(i)=sqrt(difx^2+dify^2);
theta(i)=2*(atan(dify/difx));
drawnow
delete(h)
end
theta=unwrap(theta,0.9)/2;
figure, polar(theta,rho)
np = 256; 
theta_i=linspace(0,2*pi,np);
rho_i=interp1(theta,rho,theta_i);
figure, plot(theta,rho), hold on
plot(theta_i,rho_i,'r--')


