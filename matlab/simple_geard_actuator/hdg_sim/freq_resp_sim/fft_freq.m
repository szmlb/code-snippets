main;

Ndata = 50000;

%周波数特性の描画
%実験データのグラフ化
t = time_data;
U1 = input_data';
Y1 = output_data;

%工程[0] デシメーション (200us -> 1ms)
decimate_rate = 1;
t = decimate(t, decimate_rate);
U1 = decimate(U1, decimate_rate);
Y1 = decimate(Y1, decimate_rate);

%工程[2] システム同定
dt = sampling_time * decimate_rate;

data = U1;
n=Ndata;
%n=length(data);
t=(1:n)/dt-dt;
f=t/dt/dt/n;
f=f(1:end/2+1);
fftdat_u=fft(data, n)/(n/2);
fftdat_u=fftdat_u(1:end/2+1);
%fftdat_u(1)=fftdat_u(1)/2;

data = Y1;
n=Ndata; %2のn乗を選択
%n=length(data);
t=(1:n)*dt-dt;
f=t/dt/dt/n;
f=f(1:end/2+1);
fftdat_y=fft(data, n)/(n/2);
fftdat_y=fftdat_y(1:end/2+1);
%fftdat_y(1)=fftdat_y(1)/2;

fftdat_u_mag = abs(fftdat_u);
fftdat_u_ph = angle(fftdat_u)*180/pi;

fftdat_y_mag = abs(fftdat_y);
fftdat_y_ph = angle(fftdat_y)*180/pi;

%入出力の周波数データからゲインと位相特性を描画
for i=1:1:length(fftdat_u_mag)
    fftdat_mag(i) = fftdat_y_mag(i) / fftdat_u_mag(i);
    fftdat_ph(i) = fftdat_y_ph(i) - fftdat_u_ph(i); 
end
%fftdat_ph = unwrap(fftdat_ph);

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
figure(1);
subplot(211);
%plot(f, abs(fftdat_mag)); axis tight;
semilogx(f, 20*log10(fftdat_mag)); axis tight; hold on;
semilogx(f, 20*log10(fftdat_mag_filt)); axis tight;
xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]');
xlim([10^-1,2*10^1])
grid on;
subplot(212);
%plot(f, angle(fftdat)*180/pi); axis tight;
semilogx(f, fftdat_ph); axis tight; hold on;
semilogx(f, fftdat_ph_filt); axis tight;
xlabel('Frequency [Hz]'); ylabel('Angle [degree]')
xlim([10^0,2*10^1])
grid on;
figure(2);
bode(sysdqm);