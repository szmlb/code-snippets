main;

decimate_rate = 1.0;
input_data = decimate(input_data, decimate_rate);
output_data = decimate(output_data, decimate_rate);

sampling_time = decimate_rate * sampling_time;

%matlabコマンドを用いた伝達関数特性の推定
data = iddata(output_data, input_data, sampling_time);
%data = detrend(data); %直流成分を除去

%{
%pemの初期モデルのためのノミナルモデル
Jmn = Jm/3; % / 50
Jln = Jl*2;  % Jlのロバスト安定性は高い. しかしロバスト性能は低い.
Bmn = 0.0;
Bln = 0.0;
Dsn = Ds;
Ksn = Ks*3.0; % / 3.5 % Ksの変動に対するロバスト性能は高い. しかしロバスト安定性は低い.
num_qm = s^2 + wl^2;
den_qm = Jmn * s^2 * (s^2 + wm^2);
sysqm = num_qm / den_qm;
sysdqmn = sysqm * s;
num_ql = wl^2;
den_ql = s^2 + wl^2;
sysql = num_ql / den_ql;
sysdqln = sysql * s;
%}

%pemの初期モデルのためのノミナルモデル
Jmn = Jm; % / 50
Jln = Jl;  % Jlのロバスト安定性は高い. しかしロバスト性能は低い.
Bmn = Bm;
Bln = Bl;
Dsn = Ds;
Ksn = Ks;
Ktn = 1.0;

%真の伝達関数
num_dqm = Kt * (Jl * ng^2 * s^3 + (Ds + Bl) * ng^2 * s^2 + Ks * ng^2 * s);
den_dqm = Jm * Jl * ng^2 * s^4 + ( ((Ds + Bl) * Jm + Bm * Jl) * ng^2 + Ds * Jl ) * s^3 + ((Jm * Ks + Bm * Ds + Bm * Bl) * ng^2 + Jl * Ks + Bl * Ds) * s^2 + (Bm * Ks * ng^2 + Bl * Ks) * s;
sysdqm = num_dqm / den_dqm;
sysqm = sysdqm / s;

num_dql = Kt * (Ds * ng * s + Ks * ng) * s;
den_dql = Jm * Jl * ng^2 * s^4 + ( ((Ds + Bl) * Jm + Bm * Jl) * ng^2 + Ds * Jl ) * s^3 + ((Jm * Ks + Bm * Ds + Bm * Bl) * ng^2 + Jl * Ks + Bl * Ds) * s^2 + (Bm * Ks * ng^2 + Bl * Ks) * s;

sysdql = num_dql / den_dql;
sysql = sysdql / s;

%m1 = pem(data, sysdqmn);
m1 = pem(data);
m2 = tf(m1);
m3 = d2c(m2);
%m3 = tf(m1);
m4 = spa(data); % スペクトル解析法で周波数特性を推定
%m5 = cra(data); % prewhitened-based corelation method?

freq_range_num = 1000;
w = logspace(-1, 3, freq_range_num); %[rad/s]
[m1_mag, m1_ph] = bode(m1, w);
[m2_mag, m2_ph] = bode(m2, w);
[m3_mag, m3_ph] = bode(m3, w);
[m4_mag, m4_ph] = bode(m4, w);
[nom_mag, nom_ph] = bode(sysdqm, w);
for i =1:1:freq_range_num
    
    m1_mag_(i) = m1_mag(1, 1, i);
    m2_mag_(i) = m2_mag(1, 1, i);
    m3_mag_(i) = m3_mag(1, 1, i);
    m4_mag_(i) = m4_mag(1, 1, i);
    nom_mag_(i) = nom_mag(1, 1, i);

    m1_ph_(i) = m1_ph(1, 1, i);
    m2_ph_(i) = m2_ph(1, 1, i);
    m3_ph_(i) = m3_ph(1, 1, i);
    m4_ph_(i) = m4_ph(1, 1, i);
    nom_ph_(i) = nom_ph(1, 1, i);
    
end

w_freq = w / (2.0 * pi);

figure(2);
subplot(211);
%plot(f, abs(fftdat_mag)); axis tight;
p1=semilogx(w_freq, 20*log10(m1_mag_(:, :, :))); hold on;
p1.LineWidth=2.0;
p1.LineStyle = '-.'
p2=semilogx(w_freq, 20*log10(m2_mag_(:, :, :))); hold on;
p2.LineWidth=2.0;
p2.LineStyle = '-'
p3=semilogx(w_freq, 20*log10(m3_mag_(:, :, :))); hold on;
p3.LineWidth=2.0;
p3.LineStyle = '--'
%p4=semilogx(w, 20*log10(m4_mag_(:, :, :))); hold on;
p5=semilogx(w_freq, 20*log10(nom_mag_(:, :, :))); hold on;
p5.LineWidth=2.0;
p5.LineStyle = ':'
xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]');
legend('pem','df','d2c','nom')
grid on;
subplot(212);
%plot(f, angle(fftdat)*180/pi); axis tight;
p6=semilogx(w_freq, m1_ph_); axis tight; hold on;
p6.LineStyle = '-.'
p6.LineWidth=2.0;
p7=semilogx(w_freq, m2_ph_); axis tight; hold on;
p7.LineStyle = '-'
p7.LineWidth=2.0;
p8=semilogx(w_freq, m3_ph_); axis tight; hold on;
p8.LineStyle = '--'
p8.LineWidth=2.0;
%semilogx(w, (m4_ph_)*pi/180-90); axis tight; hold on;
p9=semilogx(w_freq, nom_ph_); axis tight; hold on;
p9.LineStyle = ':'
p9.LineWidth=2.0;

xlabel('Frequency [Hz]'); ylabel('Angle [degree]')
grid on;