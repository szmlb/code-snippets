main;

Ndata = size(input_data);

%実験データのグラフ化
t = time_data;
U1 = input_data';
Y1 = output_data;

dt = sampling_time;
kt = 1.0;

%工程[0] デシメーション (200us -> 1ms)
decimate_rate = 1;
t = decimate(t, decimate_rate);
U1 = decimate(U1, decimate_rate);
Y1 = decimate(Y1, decimate_rate);

dt = dt * decimate_rate;
fs = 1/dt;

[mag, w] = tfestimate(U1, Y1, [],[],[], fs);
tfest_mag = abs(mag);
tfest_ph = angle(mag)*180/pi;

figure(1);
subplot(211);
%plot(w, abs(tfest_mag)); axis tight;
%loglog(w, tfest_mag); axis tight; hold on;
semilogx(w, 20 * log10(tfest_mag)); axis tight; hold on;
xlim([10^1, 3*10^2])
xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]');
grid on;
subplot(212);
%plot(f, angle(fftdat)*180/pi); axis tight;
semilogx(w, tfest_ph); axis tight; hold on;
xlabel('Frequency [Hz]'); ylabel('Angle [degree]')
xlim([10^1, 3*10^2])
grid on;

%figure(2);
[nom_mag, nom_ph] = bode(sysdqm, w);
for i =1:1:length(w)
    nom_mag_(i) = nom_mag(1, 1, i);
    nom_ph_(i) = nom_ph(1, 1, i);
end
subplot(211);
semilogx(w/(2.0*pi), 20 * log10(nom_mag_(:, :, :))); axis tight; hold on;
xlim([1*10^1, 2*10^2])
subplot(212);
semilogx(w/(2.0*pi), nom_ph_); axis tight; hold on;
xlim([1*10^1, 2*10^2])

%{
%ゼロ位相フィルタリング [正規化周波数でフィルタリング]
% サンプリングレートは 10kHz: 5kHzがナイキスト周波数
d = designfilt('lowpassfir','PassbandFrequency',0.05, ...
         'StopbandFrequency',0.15,'PassbandRipple',0.2, ...
         'StopbandAttenuation',65,'DesignMethod','kaiserwin');
%{
d = designfilt('lowpassiir','FilterOrder',2, ...
               'PassbandFrequency',300,'PassbandRipple',0.5, ...
               'SampleRate',1/dt)
%}          
fftdat_mag_filt = filtfilt(d,fftdat_mag);
fftdat_ph_filt = filtfilt(d,fftdat_ph);

%}