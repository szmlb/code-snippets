s = tf('s');

%% Simulation parameters definition
Jm = 0.05;
Jl = 2.5;
J = (Jm * Jl) / (Jm + Jl);
Ks = 200.0;
Ds = 0.0;
taum = 0.0;
taud = 0.0;
taul = 0.0;
sampling_time = 0.0001;

%% main loop 1[sec]
time_sample = 1 * (1 / sampling_time);
for j = 1:1:5

    %state
    xvec = [0.0; 0.0; 0.0; 0.0];
    dxvec = [0.0; 0.0; 0.0; 0.0];
    % for controller
    control_sampling_time = sampling_time;
    Jn = J;
    Jmn = Jm; % / 50
    Jln = Jl / (0.5 * j); % / 20  % Jlのロバスト安定性は高い. しかしロバスト性能は低い.
    Dsn = Ds;
    Ksn = Ks;% / (100 * j); % / 3.5 % Ksの変動に対するロバスト性能は高い. しかしロバスト安定性は低い.

    %Jlの変動にロバストなKsのノミナル値決定法
    %Ksn = Ks / Jl * Jln;
    
    %% 3DoF controller design
    %Plant transfer functions
    %P1(s)
    num_P1n = Jln * s^2 + Ksn;
    den_P1n = Jln * Jmn * s^4 + Ksn * (Jmn + Jln) * s^2;
    sysP1n = num_P1n / den_P1n;
    %P2(s)
    num_P2n = Ksn * s^2;
    den_P2n = Jln * s^2 + Ksn;
    sysP2n = num_P2n / den_P2n;
    %P3(s)
    num_P3n = 1.0;
    den_P3n = s^2;
    sysP3n = num_P3n / den_P3n;
    %Controller transfer functions
    % Gyr(s)
    tau = 0.01;%0.075
    tilde = 0.9;
    wcQ1 = 500.0;
    wcQ2 = 500.0;%800.0
    num_Gyr = tf(1.0);
    den_Gyr = (tau * s + 1.0)^4;
    sysGyr = num_Gyr / den_Gyr;
    %Q1(s)
    num_Q1 = wcQ1^2;
    den_Q1 = (s + wcQ1)^2;
    sysQ1 = num_Q1 / den_Q1;
    %Q2(s)
    num_Q2 = wcQ2^2;
    den_Q2 = (s + wcQ2)^2;
    sysQ2 = num_Q2 / den_Q2;
    %C1(s)
    sysC1 = sysGyr / (sysP1n * sysP2n * sysP3n * (1-sysQ1));
    [Ac1, Bc1, Cc1, Dc1] = ssdata(sysC1);
    %C2(s)
    sysC2 = (sysQ1 - sysQ2) / (sysP1n * (1-sysQ1));
    [Ac2, Bc2, Cc2, Dc2] = ssdata(sysC2);
    %C3(s)
    sysC3 = (sysQ2) / (sysP1n * sysP2n * (1-sysQ1));
    [Ac3, Bc3, Cc3, Dc3] = ssdata(sysC3);
    %discretization
    sysC1d = c2d(sysC1, control_sampling_time);
    [Ac1d, Bc1d, Cc1d, Dc1d] = ssdata(sysC1d);
    sysC2d = c2d(sysC2, control_sampling_time);
    [Ac2d, Bc2d, Cc2d, Dc2d] = ssdata(sysC2d);
    sysC3d = c2d(sysC3, control_sampling_time);
    [Ac3d, Bc3d, Cc3d, Dc3d] = ssdata(sysC3d);
    xvecC1 = zeros(size(Bc1), 1);
    xvecC2 = zeros(size(Bc2), 1);
    xvecC3 = zeros(size(Bc3), 1);
    
    for i = 1:1:time_sample
        time = (i - 1) * sampling_time;

    control_delay = control_sampling_time / sampling_time; %[sample]
    if mod(i, control_delay) == 0
       %% Calculation for each control sampling
       %% definition for control parameters
        if i == 1
            ;
        end
       %% feedforward and feedback controllers 
        % 3DoF controller
        if(time < 0.1)
            ql_cmd = 0.0;
        else
            ql_cmd = 1.0;
        end
        ddql_res = dxvec(4);;% + (rand - 0.5) / 1.0;
        
        %continuous controller
        %dx = Ax + Bu
        %y = Cx + Du
        % C1
        dxvecC1 = Ac1 * xvecC1 + Bc1 * ql_cmd;
        uC1 = Cc1 * xvecC1 + Dc1 * ql_cmd;    
        xvecC1 = xvecC1 + dxvecC1 * control_sampling_time;
        % C2
        dxvecC2 = Ac2 * xvecC2 + Bc2 * xvec(1); 
        uC2 = Cc2 * xvecC2 + Dc2 * xvec(1);     
        xvecC2 = xvecC2 + dxvecC2 * control_sampling_time;
        % C3
        dxvecC3 = Ac3 * xvecC3 + Bc3 * ddql_res; 
        uC3 = Cc3 * xvecC3 + Dc3 * ddql_res;     
        xvecC3 = xvecC3 + dxvecC3 * control_sampling_time;
        % control input
        taum = uC1 - uC2 - uC3;
        
        %{
        %discritized controller(うまくいってない)
        %}
        
       %% data update
        time_data(i) = time;
        xvec0_data(i) = xvec(1);
        xvec1_data(i) = xvec(2);
        xvec2_data(i) = xvec(3);
        xvec3_data(i) = xvec(4);
        cmd_data(i) = ql_cmd;
        taum_data(i) = taum;
        %""" controller end """
    end
        
    %% Plant update
    %reaction torque
    if(time < 0.3)
        taud = 0.0;
        taul = 0.0;
        %taul = 1.0 * sin(2.0 * pi * 10.0 * time);
    elseif(time < 0.4)
        taud = 0.0; %3.0
        taul = 0.0; %0.0
    else
        taud = 0.0; %3.0
        taul = 50.0; %100.0
    end
    
    % derivative calculation
    dxvec = calc_deri(xvec, Jm, Jl, Ks, Ds, taum, taud, taul);
    % euler-integration
    xvec = update(xvec, dxvec, sampling_time);

end

%% data plot
%figure(0)
plot(time_data, cmd_data, time_data, xvec0_data, time_data, xvec2_data); hold on;
end


legend('cmd', 'qm', 'ql')
grid()
xlabel('time [s]')
ylabel('ql [rad]')

