clc; clear;

% 基础参数
fs = 20e6;          % 采样率 10 MHz
T = 20e-6;          % 脉冲宽度 20 us
B = 5e6;            % 带宽 5 MHz
t = 0:1/fs:T-1/fs;  % 时间轴
N = length(t);      % 点数

% LFM信号
k = B/T;  % 斜率
x = exp(1j*pi*k*t.^2);  % 发射LFM信号（单位幅值）
x=x/sqrt(sum(x.^2));%归一化

% 匹配滤波器（时域匹配）
h = conj(fliplr(x));  % 匹配滤波器核

%接收信号的信号部分
y_s=conv(x,h);
% y_s=y_s/sqrt(sum(y_s.^2));%归一化

% SNR范围（单位dB）
SNR_dB = -20:2:20;
SNR_lin = 10.^(SNR_dB/10);

% 蒙特卡洛次数
N_MC = 1e4;

% 门限设定
P_FA_target = 0.1;  % 虚警概率目标



% 初始化统计值
P_D_sim = zeros(size(SNR_lin));
P_FA_sim = zeros(size(SNR_lin));

T_thresh_arr=zeros(1,length(SNR_lin)); % 门限值


for idx = 1:length(SNR_lin)
    snr = SNR_lin(idx);
    detect_count = 0;
    false_alarm_count = 0;
    N_sigma=sqrt(1/snr);
    % 计算门限
    T_thresh = norminv(1 - P_FA_target, 0, 1) * N_sigma; % 门限值
    T_thresh_arr(idx) = T_thresh; % 存储门限值
    for k = 1:N_MC
        % 随机生成接收信号
        % 以50%的概率生成1
        ran_num=rand(1);
        A=floor(2*ran_num);
        %信号部分
        n=(randn(1,N)+1j*randn(1,N))*N_sigma/sqrt(2);
        y_n=conv(n,h);
        y=A*y_s+y_n;

        detect_peak = max(abs(y));
        if ran_num>=0.5 && detect_peak > T_thresh
            detect_count = detect_count + 1; % 检测到目标
        elseif ran_num<0.5 && detect_peak > T_thresh
            false_alarm_count = false_alarm_count + 1; % 虚警
        end
        
    end
    % 使用一次仿真统计直接计算：
    P_D_sim(idx) = detect_count / N_MC;
    P_FA_sim(idx) = false_alarm_count / N_MC;  
end




% 理论检测概率计算
P_D_th = zeros(size(SNR_lin));
for idx = 1:length(SNR_lin)
    sqrt_snr = sqrt(SNR_lin(idx));
    T_thresh = erfcinv(2* P_FA_target) * sqrt_snr; % 门限值
    P_FA_test=0.5*(1-erf(T_thresh/sqrt_snr)); % 理论虚警概率
    P_D_th(idx) = 1-normcdf((T_thresh-1)/N_sigma,0,1); % 理论检测概率
end

f1=figure(1);
semilogy(SNR_dB, P_D_sim, 'ro-', 'DisplayName', '仿真 P_D');
hold on;
semilogy(SNR_dB, P_D_th, 'b*-', 'DisplayName', '理论 P_D');
xlabel('SNR (dB)');
ylabel('检测概率 P_D');
title(['匹配滤波器检测性能，P_{FA,target} = ', num2str(P_FA_target)]);
legend;
grid on;

P_FA_target_arr=ones(1,length(SNR_dB))*P_FA_target; % 目标虚警概率
% 虚警概率验证图
f2=figure(2);
semilogy(SNR_dB, P_FA_sim, 'k^-', 'DisplayName', '仿真 P_{FA}');
hold on;
semilogy(SNR_dB,P_FA_target_arr, 'r--', 'DisplayName', '目标 P_{FA}');
xlabel('SNR (dB)');
ylabel('虚警概率 P_{FA}');
legend;
title('虚警概率仿真统计');
grid on;
