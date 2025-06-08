
%产生线性调频信号

tau=2e-5;%线性调频信号时间长度
t_all=4e-5;%信号总时间
k=1e8;%斜率
B=k*tau;%带宽
f_c=3e4;%中心频率
fs=200*max(f_c,B);%采样率
ts=1/fs;%采样间隔
N=floor(t_all*fs)+1;%采样点数(包括端点)

% 定义采样时间点
t_sample= [-(N-1)/2:(N-1)/2]*ts;
% 从 -1 到 1，采样 2001 个点
% 卷积之后的时间点
t_new = [-(N-1):N-1] * ts ;
%频域的绘制范围。频率分辨率应该是fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % 频率范围从 -N/2 到 N/2

syms t; % 定义符号变量 t
%矩形波
A=@(t) t<tau/2&t>-tau/2;
% 线性调频信号（无相位误差）
x=@(t) A(t).*exp(1i*(pi*k*t.^2+2*pi*f_c*t));

x_sample=x(t_sample); % 采样信号

%产生回波信号

% 目标距离
R_target =1e5; % 目标距离（米）
R2_target=R_target+4e3;
%R2_target=2e5; % 目标距离（米）
c= 3e8; % 光速（米/秒）
t1= 2 * R_target / c; % 回波信号的延时（秒）
t2= 2 * R2_target / c; % 回波信号的延时（秒）

% 噪声
n1 = @(t) 0.1 * randn(size(t)); % 第一个目标的噪声
n2 = @(t) 0.1 * randn(size(t)); % 第二个目标的噪声

% 计算回波信号
s1=@(t) x(t)+n1(t);
s2=@(t) x(t+t2-t1)+n2(t+t2-t1);

% s1=@(t) x(t);
% s2=@(t) x(t+t2-t1);

s1_sample=s1(t_sample);
s2_sample=s2(t_sample);

% 通过匹配滤波器
y1=conv(s1(t_sample),fliplr(conj(x(t_sample))));
y2=conv(s2(t_sample),fliplr(conj(x(t_sample))));

y1_max= max(abs(y1)); % 第一个目标的最大值
y1_db= 20*log10(abs(y1)./y1_max); % 第一个目标的dB值
y2_db= 20*log10(abs(y2)./y1_max); % 第二个目标的dB值

% 计算对应时间坐标
t1_sample= t_new  + t1; % 第一个目标的时间坐标
t2_sample= t_new  + t2; % 第二个目标的时间坐标

% 绘制结果
f1=figure(1);
% plot(t_sample,real(x(t_sample)));
% hold on;
% 绘制3dB 线
y_3db = 20*log10(0.707); % 3dB线

% 画一条与y1_db长度相同的3dB线，x轴与t1_sample一致
plot(t1_sample, y_3db*ones(size(t1_sample)), 'k--', 'DisplayName', '3dB线');
hold on;
plot(t1_sample,abs(y1),'DisplayName', '目标1');
hold on;
plot(t2_sample,abs(y2),'DisplayName', '目标2');
title('匹配滤波器输出');
xlabel('时间 (秒)');
ylabel('幅度');
% legend('目标1','目标2');
hold off;






