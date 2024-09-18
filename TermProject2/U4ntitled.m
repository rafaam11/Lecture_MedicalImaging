xish=[0 0 0 1 0 0 0 0];
Xish=fft(xish);
figure(3)
stem(abs(Xish))