clear all, 
close all, 
clc
[audio,fs] = audioread("audiowav.wav"); 
audio=audio(44100: 110250,1);
audio=audio-mean(audio);
%%%%%%%%%%%%% Parametros para la transformada 
T = 1/fs;             % Periodo de muestreo   
L=length(audio);         % Longitud de la señal
t = (0:L-1)*T;         % Vector del tiempo
n = 2^nextpow2(L);     % Puntos de la tranformada 
f = fs*(0:(n/2))/n;    % Vector de frecuencias
%%%%%%%%%%%%% Transformada de Fourier señal uno
audioF = fft(audio,n);        % Transformada
UNO_mag = abs(audioF/L); 
UNO_angle = angle(audioF/L); 
%%%%%%%%%%%%%%%%%%
p=audioF(4).*exp(1i*2*pi*f(4).*t);
p_real=real(p);
p_imag=imag(p);
mag=sqrt(p_real.^2+p_imag.^2);
angle=(atan(p_imag./p_real))*2;
plot(angle)
angle=unwrap(angle,0.9)/2;
figure, 
ax1=subplot(1,2,1);
hold on,
xlim([-10,10])
ylim([-10,10])
ax2=subplot(1,2,2);
hold on,
ylim([-10,10])
xlim([0,t(end)])
recupera=UNO_mag(4)*L.*cos(2*pi*f(4).*t+UNO_angle(4));
for i=1:80:length(t)
subplot(1,2,1)
plot(p_imag(i),p_real(i),'.r');
[xa1,ya1] = ds2nfu(ax1,p_imag(i),p_real(i));
x=mag(i)*cos(angle(i));
y=mag(i)*sin(angle(i));
l=quiver(0,0,y,x,'b');
set(l,'AutoScaleFactor',1.01)
subplot(1,2,2)
plot(t(i),recupera(i),'.r')
[xa2,ya2] = ds2nfu(ax2,t(i),recupera(i));
xl = [xa1, xa2];
yl = [ya1,  ya2];
an=annotation('line',xl,yl);
drawnow
delete(an)
delete(l)
end





