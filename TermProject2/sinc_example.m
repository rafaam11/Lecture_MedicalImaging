
x = -20:0.01:20;
y = sinc(x);

figure(1);
plot(x,y); grid on;


fft_shift = fftshift(fft(y));
fft_shift_abs = abs(fft_shift);

figure(2);
plot(fft_shift_abs); grid on;