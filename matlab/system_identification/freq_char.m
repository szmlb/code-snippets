s = tf('s')

ng = 100;
Jm = 0.14;
Jl = 0.01;
Bm = 6.8;
Bl = 0.26;
Ks = 60.0;
Ds = 0.01;
Kt = 30.0 * 10e-3;
 
tau_i = 200 * 10e-9; %200 * 10e-9

%Plant transfer function
%i to taus
A = Jm * Jl * s^4 + ( (Ds + Bl) * Jm + Ds * Jl + Bm * Jl) * s^3 + ((Jm + Jl) * Ks + (Bm + Bl) * Ds + Bm * Bl) * s^2 + (Bm + Bl) * Ks * s;
sysP = Kt * ((Ds * s + Ks) * (Jl * s^2 + Bl * s)) / A;

%curent control loop characteristics
sysIcnt = (1 / (1 + tau_i * s));

sysReal = sysP * sysIcnt;

%% P control frequency characterisitics
Kp = 10000;
Lreal_P = Kp * sysReal %P制御時の開ループ特性
Hreal_P = (Kp * sysReal) / (1 + sysReal * Kp); %P制御時の閉ループ特性

%Frequency domain characteristics plot
figure;
subplot(221)
bode(Lreal_P)    %開ループ伝達関数のボード線図
title('開ループのボード線図')
grid on;
subplot(222)
bode(Hreal_P)    %P制御時の閉ループ伝達関数のボード線図(制御帯域の確認)
title('閉ループのボード線図')
grid on;
subplot(223)

% 正の周波数のナイキスト線図のみを表示する
% http://jp.mathworks.com/help/control/ref/nyquistoptions.html?refresh=true
Pset = nyquistoptions;
Pset.ShowFullContour = 'off'; 
nyquist(Lreal_P, Pset) %開ループ伝達関数のナイキスト線図
grid on;
subplot(224)
rlocus(sysReal)  %P制御時の根軌跡
grid on;

figure;
step(Hreal_P);   %P制御時のステップ応答
grid on;