%产生线性调频信号
c=3e8; % 光速（米/秒）
tau=2e-8;%线性调频信号时间长度
k=2e17;
B=k*tau;%带宽
lambda=3e-3;
f_c=c/lambda;%中心频率
fs=5*max(f_c,B);%采样率
ts=1/fs;%采样间隔

%矩形波
A=@(t) t<tau/2&t>-tau/2;
% 线性调频信号（无相位误差）
x=@(t) A(t).*exp(1i*(pi*k*t.^2+2*pi*f_c*t));



R1_target = 1e5; % 第一个目标距离（米）
delta_target=3;
R2_target = 1e5+delta_target; % 第二个目标距离（米）
c= 3e8; % 光速（米/秒）
tau1 = 2 * R1_target / c; % 第一个目标的回波信号延时（秒）
tau2 = 2 * R2_target / c; % 第二个目标的回波信号延时（秒）
% 噪声
n = @(t) 0.1 * randn(size(t)); % 第一个目标的噪声
% 计算混合回波信号
s = @(t) x(t-tau1) + x(t-tau2) ; % 两个目标的回波信号
%s1信号横坐标
N_edge=floor((tau2-tau1)/ts/4);
N_s_floor=floor(tau1/ts-N_edge); % 采样点数
N_s_top=floor(tau2/ts+N_edge); % 采样点数
s_t_sample=[N_s_floor:N_s_top]*ts; % 从0开始采样

t_all=2*tau;
N_x=floor(t_all*fs)+1;%采样点数(包括端点)

% 定义采样时间点
x_t_sample= [-(N_x-1)/2:(N_x-1)/2]*ts;

y_t_sample= [-(N_x-1)/2+N_s_floor:(N_x-1)/2+N_s_top]*ts;

% 通过匹配滤波器
match_filter= @(t) fliplr(conj(x(t))); % 匹配滤波器
y = conv(s(s_t_sample), match_filter(x_t_sample)); % 计算匹配滤波器输出
% 计算匹配滤波器输出的最大值
y_max = max(abs(y)); % 两个目标的最大值
y_normalized = abs(y) / y_max; % 归一化
y_db= 20*log10(y_normalized); % 转换为dB

f1=figure(1);
% 绘制匹配滤波器输出
plot(y_t_sample, abs(y), 'DisplayName', '匹配滤波器输出');
title('匹配滤波器输出 |y|');
xlabel('时间 (秒)');
ylabel('幅度');

f2=figure(2);
plot(y_t_sample, y_normalized, 'DisplayName', '匹配滤波器输出 (dB)');
xlabel('时间 (秒)');
ylabel('比值');

f3=figure(3);
plot(y_t_sample*c, y_db, 'DisplayName', '匹配滤波器输出 (dB)');
xlabel('距离（米）');
ylabel('幅度 (dB)');
title('匹配滤波器输出 (dB)');



% f4=figure(4);
% x_test=conv(x(x_t_sample),match_filter(x_t_sample)); % 计算匹配滤波器输出
% plot( abs(x_test), 'DisplayName', '匹配滤波器输出');