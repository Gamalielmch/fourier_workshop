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
figure, subplot(1,2,1)        % Gráficas
stem(f(1:10),audio_mag(1:10)), title('Magnitud') 
subplot(1,2,2)
stem(f(1:10),audio_angle(1:10)), title('Fase')


%% %%% Fasor, Real,Imaginario, Magnitud y fase
harmonic=4;
p=audioF(harmonic).*exp(-1i*2*pi*f(harmonic).*t);
p_real=real(p);
p_imag=imag(p);
[angle,mag]=cart2pol(p_real,p_imag);


%% %%% Gráficas fasores 
figure, 
hold on
axis equal
for i=1:100:length(t)
plot(p_real(i),p_imag(i),'.r')
x=mag(i)*cos(angle(i));
y=mag(i)*sin(angle(i));
l=quiver(0,0,x,y,'b');
xlim([-13,13])
ylim([-13,13])
drawnow
delete(l)
end

figure, view([30,40]),hold on, grid on, 
ylim([-20,20])
zlim([-20,20])
xlim([0,2])
for i=1:100:length(t)
plot3(t(i),p_real(i),p_imag(i),'.r')
x=mag(i)*cos(angle(i));
y=mag(i)*sin(angle(i));
l=quiver3(0,0,0,t(i),x,y,'b');
drawnow
delete(l)
end



