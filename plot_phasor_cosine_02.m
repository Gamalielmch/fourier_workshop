%% %%%%%%% Limpiar
clear all, 
close all, 
clc

%% %%%%%%% cargar señal
[audio,fs] = audioread("audiowav.wav"); 
audio=audio(44100: 110250,1);
audio=audio-mean(audio);

%% %%%%%%%%%%% Parametros para la transformada 
T = 1/fs;              % Periodo de muestreo   
L=length(audio);       % Longitud de la señal
t = (0:L-1)*T;         % Vector del tiempo
n = 2^nextpow2(L);     % Puntos de la tranformada 
f = fs*(0:(n/2))/n;    % Vector de frecuencias

%% %%%%%%%%%%% Transformada de Fourier
audioF = fft(audio,n);        % Transformada
audio_mag = abs(audioF/L);      % Magnitud
audio_angle = angle(audioF/L);  % Phase




%% %%% Fasor, Real,Imaginario, Magnitud y fase
harmonic=4;
p=audioF(harmonic).*exp(-1i*2*pi*f(harmonic).*t);
p_real=real(p);
p_imag=imag(p);
[angle,mag]=cart2pol(p_real,p_imag);

%% %%% Gráficas identidad de Euler fasores/sinusoides
figure, 
ax1=subplot(1,2,1);
hold on,
xlim([-13,13])
ylim([-13,13])
ax2=subplot(1,2,2);
hold on,
ylim([-13,13])
xlim([0,t(end)])
set(gcf ,'units','normalized')
set(gcf ,'position',[0.0702    0.5057    0.5958    0.3962])
recupera=audio_mag(harmonic)*L.*cos(2*pi*f(harmonic).*t-audio_angle(harmonic));

for i=1:150:length(t)
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





