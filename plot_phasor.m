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
figure, title('Magnitud')
stem(f(1:10),UNO_mag(1:10))
figure, title('fase')
stem(f(1:10),UNO_angle(1:10))
%%%%%%%%%%%%%%%%%%
p=audioF(4).*exp(-1i*2*pi*f(4).*t);
p_real=real(p);
p_imag=imag(p);
mag=sqrt(p_real.^2+p_imag.^2);
angle=(atan(p_imag./p_real))*2;
plot(angle)
angle=unwrap(angle,0.9)/2;
figure, hold on
xlim([-20,20])
ylim([-20,20])
axis equal
% for i=1:50:length(t)
% plot(p_real(i),p_imag(i),'.r')
% x=mag(i)*cos(angle(i));
% y=mag(i)*sin(angle(i));
% l=quiver(0,0,x,y,'b');
% drawnow
% delete(l)
% end
figure, view([50,18]),hold on, grid on, 
ylim([-20,20])
zlim([-20,20])
xlim([0,2])
for i=1:50:length(t)
plot3(t(i),p_real(i),p_imag(i),'.r')
x=mag(i)*cos(angle(i));
y=mag(i)*sin(angle(i));
l=quiver3(0,0,0,t(i),x,y,'b');
drawnow
delete(l)
end



