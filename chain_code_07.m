%% %%%%%%% Limpiar
clear all, 
close all, 
clc
%% %%%%%%% cargar imagen, binarizar y preparar imagen
I=imread('C:\Users\coyul\OneDrive\Documentos\MATLAB\im.png');
bw=imbinarize(I);
bw=imresize(bw,0.2);
bw=bwmorph(bw,'open');
bw=imfill(bw,'holes');
bw=padarray(bw,[2 2],0,'both');
bwp=bwperim(bw);
bwp=bwmorph(bwp,'shrink');


%% %%%%%%% c�digo de cadena 

[yo,xo]=find(bwp==1,1);
lon=length(find(bwp==1));
x=xo; y=yo;
neigh=zeros(8,1);
chaincode=zeros(lon,1);
for i=1:lon-1
    
    pos=[y,x+1;y-1,x+1;y-1,x;y-1,x-1;y,x-1;y+1,x-1;y+1,x;y+1,x+1];
    for j=1:8
        neigh(j)=bwp(pos(j,1),pos(j,2));
    end
    [dir,~]=find(neigh==1,1);
    chaincode(i)=dir-1;
    bwp(y,x)=0;
    y=pos(dir,1);
    x=pos(dir,2);
    imshow(bwp)
    drawnow
end

%% %%%%%%% Recuperaci�n por c�digo de cadena

x=xo; y=yo;
yc=zeros(1,lon);
xc=zeros(1,lon);
coord=[0,1;-1,1;-1,0;-1,-1;0,-1;1,-1;1,0;1,1];
for i=1:lon-1
   bwp(y,x)=1;
   yc(i)=y; xc(i)=x;
   y=y+coord(chaincode(i)+1,1);
   x=x+coord(chaincode(i)+1,2);
   imshow(bwp)
   drawnow
end
%% %%%%%%% Gr�fica cadena x e y
figure
plot(yc)
hold on 
plot(xc)