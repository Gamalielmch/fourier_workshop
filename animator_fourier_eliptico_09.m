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
cx=A0;
cy=C0;
nor=0;
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


%% Generando elipses

A0=0;
C0=0;
elip=zeros(mn+1,2,nc);
for eli=1:nc
    for nr=1:nc
        recon=zeros(mn,2);
        for t = 1 : mn
            x = 0.0;
            y = 0.0;
            for i = eli : eli
                x = x + (a(i) * cos(2 * i * pi * t / mn) + b(i) * sin(2 * i * pi * t / mn));
                y = y + (c(i) * cos(2 * i * pi * t / mn) + d(i) * sin(2 * i * pi * t / mn));
            end
            recon(t,1) = A0 + x;
            recon(t,2) = C0 + y;
        end
        recon = [recon; recon(1,1) recon(1,2)];
    end
    elip(:,:,eli)=recon;
end
%%  Gráficas elipses
fig=figure('Color',[1 1 1]);
si=get(fig,'position');
ax3=axes('Position',[0 0 1 1],'Visible','off');
set(fig,'position',[si(1)-.8*si(1),si(2)-.6*si(2),si(4)*2.5,si(4)*1.4])
sub1 =subplot(1,2,1);
pos1=get(gca,'position');
%axis off
sub2 =subplot(1,2,2);
pos2=get(gca,'position');
axes(sub2), hold on
%axis off
xlim manual
ylim manual
xlim(sub2,[-180,180])
ylim(sub2,[-160,180])
polyin = polyshape(deltaxa,deltaya);
% [cx,cy] = centroid(polyin);
plot(deltaxa-cx, deltaya-cy,'color',[0,0,0]);
pne=50;
rec=zeros(mn,2);
for n=1:mn
    rec(n,1)=sum(elip(n,1,1:pne));
    rec(n,2)=sum(elip(n,2,1:pne));
end
[xa2,ya2] = ds2nfu(sub2,rec(:,1),rec(:,2));

axes(sub1), hold on
xlim manual
ylim manual
xlim(sub1,[-180,180])
ylim(sub1,[-160,180])
hold on
color=lines(pne);
c=[0,0];
cn=pne-1;
p=cn;
oo=1;
for n=1:4:mn
    
    a=plot(sub1,elip(:,1,1), elip(:,2,1), 'color',color(1,:), 'linewidth', 2);
    b=plot(sub1,[0,elip(n,1,1)],[0, elip(n,2,1)], 'color',color(1,:), 'linewidth', 2);
    for ne=2:pne
        c(ne-1)=  plot(elip(:,1,ne)+sum(elip(n,1,1:ne-1)), elip(:,2,ne)+sum(elip(n,2,1:ne-1)), 'color',color(ne,:), 'linewidth',1);  
        p(ne-1)=plot(sub1,[sum(elip(n,1,1:ne-1)),elip(n,1,ne)+sum(elip(n,1,1:ne-1))],...
            [sum(elip(n,2,1:ne-1)),elip(n,2,ne)+sum(elip(n,2,1:ne-1)) ], 'color',color(ne,:), 'linewidth', 1);
    end
    plot(sub1, rec(n,1), rec(n,2),'.','MarkerSize',2,'color',[0,0,1])
    [xa,ya] = ds2nfu(rec(n,1),rec(n,2) );
    an=annotation('line',[xa,xa2(n)],[ya,ya2(n)]);
    
 
    plot(sub2,rec(n,1),rec(n,2),'.','MarkerSize',2,'color',[0,0,1])
    drawnow
%     if mod(n,4)==0
%         print(['C:\Users\Usuario\Documents\togif\', num2str(oo),'.png'],'-dpng') %%% Change Path
%         oo=oo+1;
%     end
    delete(p)
    delete(a)
    delete(b)
    delete(c)
    delete(an)

    
end



















% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L =length(fx);     % Length of signal
% NFFT=length(fx);
% FX = fft(fx,NFFT)/L;
% FY = fft(fy,NFFT)/L;
% magx=2*abs(FX(2:25));
% phax=angle(FX(2:25))*180/pi;
% magy=2*abs(FY(2:25));
% phay=angle(FY(2:25))*180/pi;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%5
% color=lines(length(magx)); %color harmonics
% data=360*2; %length data set
% ci=1; % Ciclos for animation
% N=linspace(0,359*ci,data); %N angles (degrees)
% oo=0;
% rn = randi([0, 1], [1, length(N)]); % Save images random rn(1)=1; rn(end)=1;
% %%%% Setup figure
% fig=figure('Color',[1 1 1]);
% si=get(fig,'position');
% ax3=axes('Position',[0 0 1 1],'Visible','off');
% set(fig,'position',[si(1)-.8*si(1),si(2)-.6*si(2),si(4)*2.5,si(4)*1.4])
% sub1 =subplot(2,2,1);
% pos1=get(gca,'position');
% axis off
% sub2 =subplot(2,2,2);
% pos2=get(gca,'position');
% sub3 =subplot(2,2,3);
% pos3=get(gca,'position');
% sub4=subplot(2,2,4);
% pos4=get(gca,'position');
% fx=fx-min(fx); fy=fy-min(fy);
% plot(fy/max(fy),fx/max(fx))
% %%%%plot circles ilustration phasor
% th = 0:pi/180:2*pi;
% magmax=sum(abs(magx));
% magmay=sum(abs(magy));
% maxmax=max(magmay);
% % xlim(gca,[0,N(end)])
% % ylim(gca,[-(magmax+magmax*.01),(magmax+magmax*.01)])
% % plot(linspace(0,359*ci,length(fx)*ci),(repmat(fx-abs(FX(1)),ci)),'r-','LineWidth',0.01)
% % hold on
% ss=0;
% jj=1;
% xo=zeros(1,length(magx)+1); yo=xo;
% xo_x(1)=0; yo_x(1)=0;
% xo_y(1)=0; yo_y(1)=0;
% xmax=[0,0];
% ymax=[0,0];
% xmin=[1,1];
% ymin=[1,1];
% for j=1:length(N)
%     for i=1:length(magx)
%         xo_x(i+1)=xo_x(i)+magx(i)*cosd(i*N(j)+phax(i)+90);
%         yo_x(i+1)=yo_x(i)+magx(i)*sind(i*N(j)+phax(i)+90);
%         xo_y(i+1)=xo_y(i)+magy(i)*cosd(i*N(j)+phay(i)+0);
%         yo_y(i+1)=yo_y(i)+magy(i)*sind(i*N(j)+phay(i)+0);
%     end
%     if xmax(1)<xo_x(i+1)
%         xmax(1)=xo_x(i+1);
%     end
%     if xmax(2)<xo_y(i+1)
%         xmax(2)=xo_y(i+1);
%     end
%     if ymax(1)<yo_x(i+1)
%         ymax(1)=yo_x(i+1);
%     end
%     if ymax(2)<yo_y(i+1)
%         ymax(2)=yo_y(i+1);
%     end
%     if xmin(1)>xo_x(i+1)
%         xmin(1)=xo_x(i+1);
%     end
%     if xmin(2)>xo_y(i+1)
%         xmin(2)=xo_y(i+1);
%     end
%     if ymin(1)>yo_x(i+1)
%         ymin(1)=yo_x(i+1);
%     end
%     if ymin(2)>yo_y(i+1)
%         ymin(2)=yo_y(i+1);
%     end
%     
%     
% end
% 
% 
% xlim(sub3,[xmin(1),xmax(1)])
% ylim(sub3,[ymin(1),ymax(1)])
% hold on
% xlim(sub2,[xmin(2),xmax(2)])
% ylim(sub2,[ymin(2),ymax(2)])
% hold on
% 
% % xlim(sub3,[-(magmax+ magmax*.000001),(magmax+magmax*.0000001)])
% % ylim(sub3,[-(magmax+ magmax*.000001),(magmax+magmax*.0000001)])
% % hold on
% % xlim(sub2,[-(magmay+magmay*.000001),(magmay+magmay*.0000001)])
% % ylim(sub2,[-(magmay+magmay*.000001),(magmay+magmay*.000001)])
% % hold on
% 
% axes(sub3), hold on
% axes(sub2), hold on
% 
% for j=1:length(N)
%     ci=1;
%     
%     for i=1:length(magx)
%         xunit = magx(i) *cos(th)+xo_x(i);
%         yunit = magx(i) *sin(th)+yo_x(i);
%         cir(i)=plot(sub3,xunit,  yunit,'color',color(ci,:));
%         lin(i)=plot(sub3,[xo_x(i),xo_x(i)+magx(i)*cosd(i*N(j)+phax(i)+90)], [yo_x(i),yo_x(i)+magx(i)*sind(i*N(j)+phax(i)+90)],'color',color(ci,:));
%         
%         
%         xunit = magy(i) *cos(th)+xo_y(i);
%         yunit = magy(i) *sin(th)+yo_y(i);
%         
%         cir2(i)=plot(sub2, xunit,  yunit,'color',color(ci,:));
%         lin2(i)=plot(sub2,[xo_y(i),xo_y(i)+magy(i)*cosd(i*N(j)+phay(i))],[yo_y(i),yo_y(i)+magy(i)*sind(i*N(j)+phay(i))],'color',color(ci,:));
%         
%         
%         xo_x(i+1)=xo_x(i)+magx(i)*cosd(i*N(j)+phax(i)+90);
%         yo_x(i+1)=yo_x(i)+magx(i)*sind(i*N(j)+phax(i)+90);
%         
%         xo_y(i+1)=xo_y(i)+magy(i)*cosd(i*N(j)+phay(i)+0);
%         yo_y(i+1)=yo_y(i)+magy(i)*sind(i*N(j)+phay(i)+0);
%         
%         if i==length(magx)
%             [xa,ya] = ds2nfu(sub3,xo_x(i+1),yo_x(i+1));
%             [xa2,ya2] = ds2nfu(sub2,xo_y(i+1),yo_y(i+1));
%             
%             
%             an=annotation('line',[xa,xa2],[ya,ya],'color',color(ci,:));
%             an2=annotation('line',[xa2,xa2],[ya,ya2],'color',color(ci,:));
%             annotation('line',[xa2,xa2+.001],[ya,ya],'color',[1 0 0]);
%             
%         end
%         if magx(i)>0
%             ci=ci+1;
%         end
%         drawnow
%     end
%     if rn(j)==1
%         print(['C:\Users\Gama\Documents\togif\', num2str(oo),'.png'],'-dpng') %%% Change Path
%         oo=oo+1;
%     end
%     delete(cir)
%     delete(lin)
%     delete(cir2)
%     delete(lin2)
%     delete(an);
%     delete(an2);
% end
% 






