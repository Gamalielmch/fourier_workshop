%Animator Phasor 2d
clear
% clc
I=imread('C:\Users\Usuario\Documents\MATLAB\phasor_animated\a1.jpg');
bw=im2bw(I(:,:,1));
%bw=fliplr(imrotate(bw,0));
nc=50;
[a,b,c,d,A0,C0,mn,dx,dy,xo,yo]=elifous(bw,nc,1);

A0=0;
C0=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dibujando elipses
mn;
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
axis off
xlim manual
ylim manual
xlim(sub2,[-300,400])
ylim(sub2,[-300,400])
polyin = polyshape(dx,dy);
[cx,cy] = centroid(polyin);
plot(dx-cx, dy-cy,'color',[0,0,0]);
pne=2;
rec=zeros(mn,2);
for n=1:mn
    rec(n,1)=sum(elip(n,1,1:pne));
    rec(n,2)=sum(elip(n,2,1:pne));
end
% plot(rec(:,1),rec(:,2),'b')
[xa2,ya2] = ds2nfu(sub2,rec(:,1),rec(:,2));

axes(sub1), hold on
xlim manual
ylim manual
xlim(sub1,[-300,400])
ylim(sub1,[-300,400])
hold on
% figure
pne=5;
color=lines(pne);
c=[0,0];
cn=pne-1;
p=cn;
oo=1;
%rec=zeros(mn,2);
for n=1:3:mn
    
    a=plot(sub1,elip(:,1,1), elip(:,2,1), 'color',color(1,:), 'linewidth', 2);
    b=plot(sub1,[0,elip(n,1,1)],[0, elip(n,2,1)], 'color',color(1,:), 'linewidth', 2);
    for ne=2:pne
        c(ne-1)=  plot(elip(:,1,ne)+sum(elip(n,1,1:ne-1)), elip(:,2,ne)+sum(elip(n,2,1:ne-1)), 'color',color(ne,:), 'linewidth',1);  
        p(ne-1)=plot(sub1,[sum(elip(n,1,1:ne-1)),elip(n,1,ne)+sum(elip(n,1,1:ne-1))],...
            [sum(elip(n,2,1:ne-1)),elip(n,2,ne)+sum(elip(n,2,1:ne-1)) ], 'color',color(ne,:), 'linewidth', 1);
    end
    plot( rec(n,1), rec(n,2),'.','MarkerSize',2,'color',[0,0,1])
    [xa,ya] = ds2nfu(rec(n,1),rec(n,2) );
    an=annotation('line',[xa,xa2(n)],[ya,ya2(n)]);
    
    axes(sub2)
    plot(rec(n,1),rec(n,2),'.','MarkerSize',2,'color',[0,0,1])
    drawnow
    if mod(n,4)==0
        print(['C:\Users\Usuario\Documents\togif\', num2str(oo),'.png'],'-dpng') %%% Change Path
        oo=oo+1;
    end
    axes(sub1)
    delete(p)
    delete(a)
    delete(b)
    delete(c)
    delete(an)
    
    %      axes(sub1), hold off
    
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






