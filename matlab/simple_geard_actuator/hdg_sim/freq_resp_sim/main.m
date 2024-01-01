s = tf('s');

%% Simulation parameters definition
%{
Jm = 0.001;
Jl = 0.01;
Bm = 0.0;
Bl = 0.0;
Ks = 100.0;
Ds = 0.0;
taum = 0.0;
taud = 0.0;
taul = 0.0;
%}
ng = 100.0;
Jm = 1.39 * 10^(-5); % / 50
%Jm = 50.0 * 10^(-5);
%Jl = 5.147 * 10^(-1);
Jl = 0.01;
Bm = 1.0 * 10^(-5);
Bl = 1.0 * 10^(-3);
Ds = 0.016;
Ks = 0.61 * 10^4; % / 3.5 % Ksの変動に対するロバスト性能は高い. しかしロバスト安定性は低い.
Kt = 1.0;

sampling_time = 0.0001;
%state
xvec = [0.0; 0.0; 0.0; 0.0];
dxvec = [0.0; 0.0; 0.0; 0.0];
% for controller
control_sampling_time = sampling_time * 1;
Jmn = Jm; % / 50
Jln = Jl;  % Jlのロバスト安定性は高い. しかしロバスト性能は低い.
Bmn = Bm;
Bln = Bl;
Dsn = Ds;
Ksn = Ks; % / 3.5 % Ksの変動に対するロバスト性能は高い. しかしロバスト安定性は低い.

resolution_m = 3000*4;
resolution_l = 100000;
resolution_max_num = 2.0 * pi;

time = 10.0;
time_sample = time * (1 / sampling_time);

%真の伝達関数
num_dqm = Kt * (Jl * ng^2 * s^3 + (Ds + Bl) * ng^2 * s^2 + Ks * ng^2 * s);
den_dqm = Jm * Jl * ng^2 * s^4 + ( ((Ds + Bl) * Jm + Bm * Jl) * ng^2 + Ds * Jl ) * s^3 + ((Jm * Ks + Bm * Ds + Bm * Bl) * ng^2 + Jl * Ks + Bl * Ds) * s^2 + (Bm * Ks * ng^2 + Bl * Ks) * s;
sysdqm = num_dqm / den_dqm;
sysqm = sysdqm / s;

num_dql = Kt * (Ds * ng * s + Ks * ng) * s;
den_dql = Jm * Jl * ng^2 * s^4 + ( ((Ds + Bl) * Jm + Bm * Jl) * ng^2 + Ds * Jl ) * s^3 + ((Jm * Ks + Bm * Ds + Bm * Bl) * ng^2 + Jl * Ks + Bl * Ds) * s^2 + (Bm * Ks * ng^2 + Bl * Ks) * s;

sysdql = num_dql / den_dql;
sysql = sysdql / s;

%同定用M系列信号
input_data = idinput(20000, 'prbs', [0, 1],[-4, 4]);

%% main loop 1[sec]
for i = 1:1:time_sample
    time = (i - 1) * sampling_time;
    control_delay = control_sampling_time / sampling_time; %[sample]
    if (mod(i, control_delay) == 0 || i == 1 )
       %% Calculation for each control sampling
       %% definition for control parameters
        if i == 1
            qm_res = 0.0;
            ql_res = 0.0;
            dqm_res = 0.0;
            dql_res = 0.0;
        end
        
        qm_res = quantization(xvec(1), resolution_m, resolution_max_num); 
        dqm_res = xvec(2);
        ql_res = quantization(xvec(3), resolution_l, resolution_max_num);        
        dql_res = xvec(4);
        
        %{
        %M系列信号を用いる場合
        if i <= length(input_data)
            taum = input_data(i);
        elseif i < 2 * length(input_data)
            taum = input_data(i - length(input_data));
        else
            taum = input_data(i - length(input_data) * 2);
        end
        %}
        
        %sweep関数を用いる場合
        freq_max = 100.0;
        amp_max = 0.5;
        freq_charp =  freq_max / (time_sample * sampling_time) * time;
        taum = (0.1 + amp_max / (time_sample * sampling_time) * time) * sin(2.0 * pi * freq_charp * time);  
        
        %""" controller end """
    end
    
    %% data update
    time_data(i) = time;
    xvec0_data(i) = xvec(1);
    xvec1_data(i) = xvec(2);
    xvec2_data(i) = xvec(3);
    xvec3_data(i) = xvec(4);
    dxvec3_data(i) = dxvec(4);
    qm_data(i) =  qm_res;
    dqm_data(i) =  dqm_res;
    ql_data(i) =  ql_res;
    dql_data(i) =  dql_res;
    taum_data(i) = taum;
        
    input_data(i) = taum;
    output_data(i) = dqm_res;
        
    %% Plant update
    %reaction torque
    if(time < 0.3)
        taud = 0.0;
        taul = 0.0;
        %taul = 1.0 * sin(2.0 * pi * 10.0 * time);
    elseif(time < 0.3)
        taud = 0.0; %3.0
        taul = 0.0; %0.0
    else
        taud = 0.0; %3.0
        taul = 0.0;%50.0; %100.0
    end
    
    % derivative calculation
    %type = 0; %type1(ギアなし)
    %type = 1; %type 2(ギア付き)    
    %dxvec = calc_deri(type, xvec, Jm, Jl, Bm, Bl, Ks, Ds, ng, taum, taud, taul);
    %xvec = update(xvec, dxvec, sampling_time);
    % runge-kutta integration
    %type = 0; %type1
    type = 1; %type 2
    param = [Jm, Jl, Bm, Bl, Ks, Ds, ng, taum, taud, taul];
    xvec = rk4_update(type, xvec, time, sampling_time, param);
   
end

%中央差分
dt = sampling_time;
for i=1:1:length(qm_data)
    if i == 1 dqm_data(i) = (qm_data(i+1)-qm_data(i)) / dt;
    elseif i == length(qm_data) dqm_data(i) = (qm_data(i)-qm_data(i-1)/dt);
    else dqm_data(i) = (qm_data(i+1)-qm_data(i-1)) / (2.0*dt);
    end
end
dqm_data(length(qm_data)) = dqm_data(length(qm_data)-1);
for i=1:1:length(ql_data)
    if i == 1 dql_data(i) = (ql_data(i+1)-ql_data(i)) / dt;
    elseif i == length(ql_data) dql_data(i) = (ql_data(i)-ql_data(i-1)/dt);
    else dql_data(i) = (ql_data(i+1)-ql_data(i-1)) / (2.0*dt);
    end
end
dql_data(length(ql_data)) = dql_data(length(ql_data)-1);

%同定用の処理
output_data = dqm_data;
%input_data = input_data'; sinスイープの時は使用する
output_data = output_data';

%% data plot
figure(1)
subplot(4, 1, 1)
plot(time_data, xvec0_data, time_data, qm_data)
legend('real', 'quant')
grid()
xlabel('time [s]')
ylabel('angle [rad]')
subplot(4, 1, 2)
plot(time_data, xvec1_data, time_data, dqm_data)
legend('real', 'quant')
grid()
xlabel('time [s]')
ylabel('angular velocity [rad/s]')
figure(1)
subplot(4, 1, 3)
plot(time_data, xvec2_data, time_data, ql_data)
legend('real', 'quant')
grid()
xlabel('time [s]')
ylabel('angle [rad]')
subplot(4, 1, 4)
plot(time_data, xvec3_data, time_data, dql_data)
legend('real', 'quant')
grid()
xlabel('time [s]')
ylabel('angular velocity [rad/s]')

% data plot
input_data = input_data(1:time_sample);
figure(2)
plot(time_data, input_data');
grid()
xlabel('time [s]')
ylabel('torque [Nm]')

%data write
delimiterOut = ' ';
dlmwrite('./time.dat', time_data', delimiterOut);
dlmwrite('./input.dat', taum_data', delimiterOut);
dlmwrite('./qm.dat', qm_data', delimiterOut);
dlmwrite('./ql.dat', ql_data', delimiterOut);