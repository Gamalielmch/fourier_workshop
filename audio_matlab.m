clear all, 
close all, 
clc
[audio,fs] = audioread("audiowav.wav"); 
%sound(audio,fs)
audio=audio(44100:44900,1);
%%%%%%%%%%%%% Parametros para la transformada 
T = 1/fs;             % Periodo de muestreo   
L=length(audio);         % Longitud de la señal
t = (0:L-1)*T;         % Vector del tiempo
n = 2^nextpow2(L);     % Puntos de la tranformada 
f = fs*(0:(n/2))/n;    % Vector de frecuencias
% audio=50*cos(2*pi*1000.*t);
%%%%%%%%%%%%% Transformada de Fourier señal uno
audioF = fft(audio,n);        % Transformada
UNO_mag = abs(audioF/L); 
UNO_angle = angle(audioF/L); 
figure 
plot(f,UNO_mag(1:(n/2)+1))
figure
plot(f,UNO_angle(1:(n/2)+1))
xr=audio'*0;
figure, stem(f,UNO_mag(1:513))
for i=1:20
  xr=xr+UNO_mag(i).*pi/2*cos(2*pi*f(i).*t+UNO_angle(i));
end 
figure, plot(audio), hold on 
plot(xr,'r--')