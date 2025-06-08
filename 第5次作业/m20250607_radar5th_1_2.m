
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
R_list=[0,1e3,1.5e3, 2e3,2.5e3,3e3, 4e3, 8e3]; % 目标距离列表
%R2_target=2e5; % 目标距离（米）
c= 3e8; % 光速（米/秒）
tau= 2 * R_target / c; % 回波信号的延时（秒）

% 噪声
n = @(t) 0.1 * randn(size(t)); % 第一个目标的噪声

% 计算回波信号
s=@(t) x(t)+n(t);

% 通过匹配滤波器
y=conv(s(t_sample),fliplr(conj(x(t_sample))));

y_max= max(abs(y)); % 第一个目标的最大值
% 3dB 线
y_3db = 20*log10(0.707); % 3dB线
% 绘制结果
f1=figure(1);

title('匹配滤波器输出');
xlabel('时间 (秒)');
ylabel('幅度');
%绘制3dB线
plot(t_sample+tau, y_3db*ones(size(t_sample)), 'k--', 'DisplayName', '3dB线');
hold on;


for index=1:1:length(R_list)
    R1_target =R_list(index);
    t1= 2 * R1_target / c;%相对时间差
    s1 =@(t) x(t+t1)+n(t+t1);
    y1=conv(s1(t_sample),fliplr(conj(x(t_sample))));
    y1_db=20*log10(abs(y1)./y_max); % 目标的dB值
    t1_sample= t_new  + t1 + tau; % 目标的时间坐标
    plot(t1_sample,y1_db,'DisplayName', ['距离' ,num2str(R1_target)]);
    hold on;
end
legend('show');
hold off;