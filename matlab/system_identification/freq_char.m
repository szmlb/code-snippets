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
Lreal_P = Kp * sysReal %P���䎞�̊J���[�v����
Hreal_P = (Kp * sysReal) / (1 + sysReal * Kp); %P���䎞�̕��[�v����

%Frequency domain characteristics plot
figure;
subplot(221)
bode(Lreal_P)    %�J���[�v�`�B�֐��̃{�[�h���}
title('�J���[�v�̃{�[�h���}')
grid on;
subplot(222)
bode(Hreal_P)    %P���䎞�̕��[�v�`�B�֐��̃{�[�h���}(����ш�̊m�F)
title('���[�v�̃{�[�h���}')
grid on;
subplot(223)

% ���̎��g���̃i�C�L�X�g���}�݂̂�\������
% http://jp.mathworks.com/help/control/ref/nyquistoptions.html?refresh=true
Pset = nyquistoptions;
Pset.ShowFullContour = 'off'; 
nyquist(Lreal_P, Pset) %�J���[�v�`�B�֐��̃i�C�L�X�g���}
grid on;
subplot(224)
rlocus(sysReal)  %P���䎞�̍��O��
grid on;

figure;
step(Hreal_P);   %P���䎞�̃X�e�b�v����
grid on;