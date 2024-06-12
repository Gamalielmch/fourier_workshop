%function [a,b,c,d,A0,C0,mn,deltaxa,deltaya,xo,yo]=elifous(BW,nc,nor)
%%función basada en Fourier elíptico desarrollado por Kuhl & Giardina
%
% varibles de entrada:
%  BW: imagen binaria del contorno
%
%  varibles de salida:
%  a:  coeficientes de los armonicos correspondiente a la proyección en x
%      del código de cadena
%  b:  coeficientes de los armonicos correspondiente a la proyección en x
%      del código de cadena
%  c:  coeficientes de los armonicos correspondiente a la proyección en y
%      del código de cadena
%  d:  coeficientes de los armonicos correspondiente a la proyección en y
%      del código de cadena
%  A0: componente DC correspondiente a la proyección en x
%  C0: componente DC correspondiente a la proyección en y
%  mn: longitud del código de cadena

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
bw=bwperim(bw);
bw=bwmorph(bw,'shrink');


bw = padarray(bw,2);
bw=bwperim(bw);
bw=bwmorph(bw,'spur','inf');
bw=bwmorph(bw,'thin','inf');
bw=bwmorph(bw,'spur','inf');
pxl= regionprops(bw,'PixelList');
pxl=pxl.PixelList;

%% %%%%%%% código de cadena 

[yo,xo]=find(bw==1,1);
lon=length(find(bw==1));
x=xo; y=yo;
neigh=zeros(8,1);
chaincode=zeros(lon,1);
%%                               Código de cadena
[yi,xi]=find(bw,1,'first');
xo=xi; yo=yi;
con=2;
longitud=length(pxl);
bw=double(bw);
le=sqrt(2);

cadena=zeros(1,longitud);
for i=1:lon-1
    pos=[y,x+1;y-1,x+1;y-1,x;y-1,x-1;y,x-1;y+1,x-1;y+1,x;y+1,x+1];
    for j=1:8
        neigh(j)=bw(pos(j,1),pos(j,2));
    end
    [dir,~]=find(neigh==1,1);
    chaincode(i)=dir-1;
    bw(y,x)=0;
    y=pos(dir,1);
    x=pos(dir,2);
end

mn=length(chaincode);


%% Tiempo, distancia Traversal y proyección de cada en eje x e y
deltat=zeros(1,length(chaincode));
deltat(1)=1+((sqrt(2)-1)/2)*(1-(-1)^chaincode(1));
deltat2=1+((sqrt(2)-1)/2)*(1-(-1).^chaincode);
deltaxa=sign(6-chaincode(1))*sign(2-chaincode(1));
deltaya=sign(4-chaincode(1))*sign(chaincode(1));
for i=2:length(chaincode)
    deltat(i)=deltat(i-1)+1+((sqrt(2)-1)/2)*(1-(-1)^chaincode(i));
    deltaxa(i)=deltaxa(i-1)+ sign(6-chaincode(i))*sign(2-chaincode(i));
    deltaya(i)=deltaya(i-1)+ sign(4-chaincode(i))*sign(chaincode(i));
end
T = deltat(end);
deltax=sign(6-chaincode).*sign(2-chaincode);
deltay=sign(4-chaincode).*sign(chaincode);

%% Cálculo de los coeficientes
nc=256;
a=zeros(1,nc);
b=a;
c=a;
d=a;

q_x = deltax./ deltat2;
q_y = deltay./ deltat2;
for n=1:nc
    two_n_pi = 2 * n * pi;
    sigma_a = 0;
    sigma_b = 0;
    sigma_c = 0;
    sigma_d = 0;
    tp_prev = 0;
    for p=1:length(cadena)
        sigma_a = sigma_a + q_x(p) * (cos(two_n_pi * deltat(p) / T) - cos(two_n_pi * tp_prev / T));
        sigma_b = sigma_b + q_x(p) * (sin(two_n_pi * deltat(p) / T) - sin(two_n_pi * tp_prev / T));
        sigma_c = sigma_c + q_y(p) * (cos(two_n_pi * deltat(p) / T) - cos(two_n_pi * tp_prev / T));
        sigma_d = sigma_d + q_y(p) * (sin(two_n_pi * deltat(p) / T) - sin(two_n_pi * tp_prev / T));
        tp_prev = deltat(p);
    end
    r = T/(2*n^2*pi^2);
    
    a(n) = r * sigma_a;
    b(n) = r * sigma_b;
    c(n)= r * sigma_c;
    d(n)= r * sigma_d;
end


%%  Cálculo de la componente directa
A0 = 0;
C0 = 0;
A0 = A0 + deltax(1) / (2 * deltat2(1)) * (deltat(1))^2 ;
C0 = C0 + deltay(1) / (2 * deltat2(1)) * (deltat(1))^2 ;
for p = 2 : length(cadena)
    
    zeta = deltaxa(p - 1) - deltax(p) / deltat2(p) * deltat(p - 1);
    delta = deltaya(p - 1) - deltay(p) / deltat2(p) * deltat(p - 1);
    A0 = A0 + deltax(p) / (2 * deltat2(p)) * ((deltat(p))^2 - (deltat(p - 1))^2) + zeta * (deltat(p) - deltat(p-1));
    C0 = C0 + deltay(p) / (2 * deltat2(p)) * ((deltat(p))^2 - (deltat(p - 1))^2) + delta * (deltat(p) - deltat(p-1));
    
end
A0 = A0 / T;
C0 = C0 / T;

nor=1;
%%  Normalization
if nor == 1
    % Remove DC components
    A0 = 0;
    C0 = 0;
    
    % Compute theta1
    theta1 = 0.5 * atan2(2 * (a(1) * b(1) + c(1) * d(1)) , ...
        (a(1)^2 + c(1)^2 - b(1)^2 - d(1)^2));
    if theta1<0
        theta1=theta1+2*pi;
    end
    
    costh1 = cos(theta1);
    sinth1 = sin(theta1);
    
    a_star_1 = costh1 * a(1) + sinth1 * b(1);
    %     b_star_1 = -sinth1 * a(1) + costh1 * b(1);
    c_star_1 = costh1 * c(1) + sinth1 * d(1);
    %     d_star_1 = -sinth1 * c(1) + costh1 * d(1);
    
    % Compute psi1
    psi1 = atan(c_star_1 / a_star_1) ;
    
    % Compute E
    E = sqrt(a_star_1^2 + c_star_1^2);
    
    cospsi1 = cos(psi1);
    sinpsi1 = sin(psi1);
    
    for i = 1 : nc
        normalized = [cospsi1 sinpsi1; -sinpsi1 cospsi1] * [a(i) b(i); c(i) d(i)] * ...
            [cos(theta1 * i) -sin(theta1 * i); sin(theta1 * i) cos(theta1 * i)];
        a(i) = normalized(1,1) / E;
        b(i) = normalized(1,2) / E;
        c(i) = normalized(2,1) / E;
        d(i) = normalized(2,2) / E;
    end
    
end
