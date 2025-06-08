PRF=1e4; % 脉冲重复频率（Hz）

%产生线性调频信号
c=3e8; % 光速（米/秒）
tau=2e-8;%线性调频信号时间长度
t_all=4e-8;%信号总时间
k=1e11;%斜率
B=k*tau;%带宽
lambda=3e-2;
f_c=c/lambda;%中心频率
% f_c=3e4;%中心频率
fs=20*max(f_c,B);%采样率
ts=1/fs;%采样间隔
N=floor(t_all*fs)+1;%采样点数(包括端点)

% 定义采样时间点
t_sample= [-(N-1)/2:(N-1)/2]*ts;
% 从 -1 到 1，采样 2001 个点
% 卷积之后的时间点
t_new = [-(N-1):N-1] * ts ;
%频域的绘制范围。频率分辨率应该是fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % 频率范围从 -N/2 到 N/2

%目标
v=30; % 目标速度（米/秒）
f_d=2*v/lambda; % 多普勒频移

syms t; % 定义符号变量 t
% 矩形窗
A = @(t) (t > -tau/2) & (t < tau/2);

% 发射信号（参考）
x_tx = @(t) A(t).*exp(1i*(pi*k*t.^2 + 2*pi*f_c*t));

% 单个回波：多普勒频偏 + 延迟
x_rx = @(t, n) A(t).*exp(1i*(pi*k*t.^2 + 2*pi*(f_c + f_d)*t))+0.1 * randn(size(t));

% 匹配滤波器模板（参考信号共轭时间反转）
ref = conj(fliplr(x_tx(t_sample)));
f1=figure(1);
x_sample = x_tx(t_sample); % 采样信号
% plot(t_sample, abs(x_sample), 'DisplayName', '实部');
% hold on;
% title('线性调频信号');
% xlabel('时间 (秒)');
% ylabel('幅度');
% legend('show');


t_accumulate = 0.1;             % 积累时间（秒）
num_pulses = floor(PRF * t_accumulate);  % 脉冲个数
compressed_pulses = zeros(num_pulses, 2*N - 1); % 储存压缩后的脉冲
for n = 0:num_pulses-1
    % 生成第n个回波
    echo = x_rx(t_sample);
    % 匹配滤波器（脉冲压缩）
    compressed = conv(echo, ref);
    compressed_pulses(n+1, :) = compressed;
end
% 从所有压缩后的波形中选择一个时间点（通常选择主峰对应的点）
[~, peak_idx] = max(mean(abs(compressed_pulses), 1));
doppler_data = compressed_pulses(:, peak_idx);  % 沿快时间做慢时间FFT

% 多普勒FFT
N_doppler = 2^nextpow2(num_pulses); % 2的幂长度
doppler_spectrum = fftshift(fft(doppler_data, N_doppler));

% 多普勒频率轴
f_d_axis = linspace(-PRF/2, PRF/2, N_doppler);
velocity_axis = f_d_axis * lambda / 2; % 由多普勒频率转为速度

% 匹配滤波二维图像
figure;
imagesc((1:num_pulses)/PRF, t_new*1e6, 20*log10(abs(compressed_pulses.')));
xlabel('时间 (秒)');
ylabel('延迟时间 (μs)');
title('匹配滤波结果（脉冲压缩输出）');
colorbar;

% 多普勒谱图
figure;
plot(velocity_axis, 20*log10(abs(doppler_spectrum)/max(abs(doppler_spectrum))));
xlabel('速度 (m/s)');
ylabel('幅度 (dB)');
title('多普勒谱（速度估计）');
grid on;


% % 脉冲压缩v
% y=conv(x(t_sample),fliplr(conj(x(t_sample))));
% y_max= max(abs(y)); % 第一个目标的最大值
% % 3dB 线
% y_3db = 20*log10(0.707); % 3dB线
% % 绘制结果
% f2=figure(2);
% title('匹配滤波器输出');
% xlabel('时间 (秒)');
% ylabel('幅度');
% plot(t_sample,y_3db*ones(size(t_sample)), 'k--', 'DisplayName', '3dB线');
% hold on;
% % 绘制匹配滤波器输出
% plot(t_new, 20*log10(abs(y)./y_max), 'DisplayName', '匹配滤波器输出');
% legend('show');
% hold off;
% % 绘制频谱
% f3=figure(3);
% title('匹配滤波器输出的频谱');
% xlabel('频率 (Hz)');
% ylabel('幅度谱');
% % 计算频谱
% Y=fftshift(fft(y,N)); % 频谱
% plot(f,abs(Y));
% legend('show');
% hold off;