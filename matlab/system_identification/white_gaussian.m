u = idinput(250001, 'rgs', [0.01, 0.99],[-1, 1]);
dlmwrite('white_gaussian.dat', u, ' ');

figure;
subplot(221), stairs(u), axis([0 100 -3.5 3.5])
subplot(222), plot(covf(u, 1000)), axis([0 1000 0 0.5])
subplot(223), periodogram(u)
subplot(224), hist(u,100)

figure;
plot(u)