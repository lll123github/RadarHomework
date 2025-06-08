%产生线性调频信号
c=3e8; % 光速（米/秒）
t_accumulate=0.1; % 累积时间
k=1e11;%斜率
B=k*t_accumulate;%带宽
lambda=3e-2;
f_c=c/lambda;%中心频率
% f_c=3e4;%中心频率
fs=20*max(f_c,B);%采样率
ts=1/fs;%采样间隔
N=floor(t_accumulate*fs)+1;%采样点数(包括端点)

% 定义采样时间点
t_sample= [-(N-1)/2:(N-1)/2]*ts;

% 请求的 20000000001x1 (149.0GB)数组超过预设的最大数组大小。创建大于此限制的数组可能需要较长时间，并且会导致 MATLAB 无响应。有关详细信息，请参阅数组大小限制或预设项面板。

% 出错 m20250608_radar5th_2_1 (line 14)
t_sample= [-(N-1)/2:(N-1)/2]*ts;


% 从 -1 到 1，采样 2001 个点
% 卷积之后的时间点
t_new = [-(N-1):N-1] * ts ;
%频域的绘制范围。频率分辨率应该是fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % 频率范围从 -N/2 到 N/2


PRF=1e4;
t_accumulate=0.1; % 脉冲串累积时间
duty_ratio=0.1; % 占空比
syms t; % 定义符号变量 t

num=PRF*t_accumulate; % 累积次数
tau=t_accumulate*duty_ratio/num; % 每次脉冲的时间长度
% 脉冲间隔时长
tau_interval=t_accumulate/PRF; % 脉冲间隔时长
%矩形波
A=@(t) t<tau & t>=0;
for i=-num/2:num/2-1
    lower= i*tau;
    upper= (i+1)*tau;
    A(t)=A(t) | (t>i*tau_interval & t<i*tau_interval+tau);
end


% 线性调频脉冲串信号（无相位误差）
x=@(t) A(t).*exp(1i*(pi*k*t.^2+2*pi*f_c*t));

x_sample=x(t_sample); % 采样信号

