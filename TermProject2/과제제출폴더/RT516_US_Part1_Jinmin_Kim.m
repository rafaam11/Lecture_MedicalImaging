%------------------------------------------------------
%   2021 RT516 - Medical Imaging
%   Term Project 2, Part 1 - Ultrasound Imaging
%   Part I : Design of Linear Array Transducers
%   202123008 KIM Jinmin (M.S. Candidate)
%   Department of Robotics Engineering
%------------------------------------------------------
clear all
close all
clc

%% Parameters
% 길이 단위는 mm로 통일
ele = 64;
len = 1000;
d = 0.2;             % mm (m * e-03)
z = 30;              % mm (m * e-03)
Freq = 3;            % MHz
lambda = 1/Freq;
D = len*d;
W = 1;
X = -D/2+d:d:D/2;     


%% Question 0 : No weighting / apodization

Y1 = zeros(1,len);           
Y1(len/2-ele/2-1:len/2+ele/2) = W;    

% Figure 1
figure;     
subplot(2,1,1);
plot(X,Y1); 
title('Aperture1 : No weighting / apodization');
axis([-D/2,D/2,0,1]); 
xlabel('mm'); ylabel('Weighting [a.u.]');

subplot(2,1,2);
fft_shift1 = fftshift(fft(Y1));
fft_shift_abs1 = abs(fft_shift1);
ydb1 = mag2db(fft_shift_abs1);
plot(X,ydb1); grid on;
title('2-way far field beam pattern 1 (PSF)');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
hold on;
yline(-6,'-r');
dB_max1 = max(ydb1);
yline(dB_max1,'-g');
hold off;
legend('PSF (dB)','-6dB line','y max (dB)');

%% Question 1 : Gaussian apodization function

L = 100;
w = gausswin(L);
Y2 = zeros(1,len);           
Y2(len/2-ele/2-1:len/2+ele/2) = w(L/2-ele/2-1:L/2+ele/2);    

% Figure 2
figure;     
subplot(2,1,1);
plot(X,Y2); 
title('Aperture2 : Gaussian apodization');
axis([-D/2,D/2,0,1]); 
xlabel('mm'); ylabel('Weighting [a.u.]');

subplot(2,1,2);
fft_shift2 = fftshift(fft(Y2));
fft_shift_abs2 = abs(fft_shift2);
ydb2 = mag2db(fft_shift_abs2);
plot(X,ydb2); grid on;
title('2-way far field beam pattern 2 (PSF)');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
hold on;
yline(-6,'-r');
dB_max2 = max(ydb2);
yline(dB_max2,'-g');
hold off;
legend('PSF (dB)','-6dB line','y max (dB)');


%% Question 2 : Spread the elements
len = 1000;
d = 0.2*2;             % mm (m * e-03)
D = len*d;
X = -D/2+d:d:D/2;     

Y3 = zeros(1,len);           
Y3(len/2-ele/2-1:len/2+ele/2) = W;    

% Figure 3
figure;     
subplot(2,1,1);
plot(X,Y3); 
title('Aperture3 : Spread the elements');
axis([-D/4,D/4,0,1]); 
xlabel('mm'); ylabel('Weighting [a.u.]');

subplot(2,1,2);
fft_shift3 = fftshift(fft(Y3));
fft_shift_abs3 = abs(fft_shift3);
ydb3 = mag2db(fft_shift_abs3);
plot(X,ydb3); grid on;
title('2-way far field beam pattern 3 (PSF)');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
hold on;
yline(-6,'-r');
dB_max3 = max(ydb3);
yline(dB_max3,'-g');
hold off;
legend('PSF (dB)','-6dB line','y max (dB)');

%% Question 3 : Design my own aperture

len = 1000;
d = 0.2;             % mm (m * e-03)
D = len*d;
X = -D/2+d:d:D/2;     

L = 100;
w = gausswin(L,1.2);
Y4 = zeros(1,len);           
Y4(len/2-ele/2-1:len/2+ele/2) = w(L/2-ele/2-1:L/2+ele/2);    

% Figure 4
figure;     
subplot(2,1,1);
plot(X,Y4); 
title('Aperture4 : Design my own aperture');
axis([-D/2,D/2,0,1]); 
xlabel('mm'); ylabel('Weighting [a.u.]');

subplot(2,1,2);
fft_shift4 = fftshift(fft(Y4));
fft_shift_abs4 = abs(fft_shift4);
ydb4 = mag2db(fft_shift_abs4);
plot(X,ydb4); grid on;
title('2-way far field beam pattern 4 (PSF)');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
hold on;
yline(-6,'-r');
dB_max4 = max(ydb4);
yline(dB_max4,'-g');
hold off;
legend('PSF (dB)','-6dB line','y max (dB)');

%% Comparison

figure;

subplot(3,1,1); 
% figure 1
len = 1000; d = 0.2; D = len*d; X = -D/2+d:d:D/2;
Y1 = zeros(1,len);           
Y1(len/2-ele/2-1:len/2+ele/2) = W;    
fft_shift1 = fftshift(fft(Y1));
fft_shift_abs1 = abs(fft_shift1);
ydb1 = mag2db(fft_shift_abs1);
plot(X,ydb1); grid on; 
hold on;
% figure 2
len = 1000; d = 0.2; D = len*d; X = -D/2+d:d:D/2;   
L = 100;
w = gausswin(L);
Y2 = zeros(1,len);           
Y2(len/2-ele/2-1:len/2+ele/2) = w(L/2-ele/2-1:L/2+ele/2);    
fft_shift2 = fftshift(fft(Y2));
fft_shift_abs2 = abs(fft_shift2);
ydb2 = mag2db(fft_shift_abs2);
plot(X,ydb2);
hold off;
% comparison
title('Figure 1 vs Figure 2');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
legend('Figure 1','Figure 2');

subplot(3,1,2);
% figure 1
len = 1000; d = 0.2; D = len*d; X = -D/2+d:d:D/2;
Y1 = zeros(1,len);           
Y1(len/2-ele/2-1:len/2+ele/2) = W;    
fft_shift1 = fftshift(fft(Y1));
fft_shift_abs1 = abs(fft_shift1);
ydb1 = mag2db(fft_shift_abs1);
plot(X,ydb1); grid on; 
hold on;
% figure 3
len = 1000; d = 0.2*2; D = len*d; X = -D/2+d:d:D/2;   
Y3 = zeros(1,len);           
Y3(len/2-ele/2-1:len/2+ele/2) = W; 
fft_shift3 = fftshift(fft(Y3));
fft_shift_abs3 = abs(fft_shift3);
ydb3 = mag2db(fft_shift_abs3);
plot(X,ydb3);
hold off;
% comparison
title('Figure 1 vs Figure 3');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
legend('Figure 1','Figure 3');


subplot(3,1,3);
% figure 1
len = 1000; d = 0.2; D = len*d; X = -D/2+d:d:D/2;
Y1 = zeros(1,len);           
Y1(len/2-ele/2-1:len/2+ele/2) = W;    
fft_shift1 = fftshift(fft(Y1));
fft_shift_abs1 = abs(fft_shift1);
ydb1 = mag2db(fft_shift_abs1);
plot(X,ydb1); grid on; 
hold on;
% figure 4
len = 1000; d = 0.2; D = len*d; X = -D/2+d:d:D/2;   
L = 100;
w = gausswin(L,1.2);
Y4 = zeros(1,len);           
Y4(len/2-ele/2-1:len/2+ele/2) = w(L/2-ele/2-1:L/2+ele/2);       
fft_shift4 = fftshift(fft(Y4));
fft_shift_abs4 = abs(fft_shift4);
ydb4 = mag2db(fft_shift_abs4);
plot(X,ydb4); 
hold off;
% comparison
title('Figure 1 vs Figure 4');
axis([-100,100,-50,50]); 
xlabel('mm'); ylabel('dB');
legend('Figure 1','Figure 4');


