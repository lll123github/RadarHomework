%�������Ե�Ƶ�ź�
c=3e8; % ���٣���/�룩
t_accumulate=0.1; % �ۻ�ʱ��
k=1e11;%б��
B=k*t_accumulate;%����
lambda=3e-2;
f_c=c/lambda;%����Ƶ��
% f_c=3e4;%����Ƶ��
fs=20*max(f_c,B);%������
ts=1/fs;%�������
N=floor(t_accumulate*fs)+1;%��������(�����˵�)

% �������ʱ���
t_sample= [-(N-1)/2:(N-1)/2]*ts;

% ����� 20000000001x1 (149.0GB)���鳬��Ԥ�����������С���������ڴ����Ƶ����������Ҫ�ϳ�ʱ�䣬���һᵼ�� MATLAB ����Ӧ���й���ϸ��Ϣ������������С���ƻ�Ԥ������塣

% ���� m20250608_radar5th_2_1 (line 14)
t_sample= [-(N-1)/2:(N-1)/2]*ts;


% �� -1 �� 1������ 2001 ����
% ���֮���ʱ���
t_new = [-(N-1):N-1] * ts ;
%Ƶ��Ļ��Ʒ�Χ��Ƶ�ʷֱ���Ӧ����fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % Ƶ�ʷ�Χ�� -N/2 �� N/2


PRF=1e4;
t_accumulate=0.1; % ���崮�ۻ�ʱ��
duty_ratio=0.1; % ռ�ձ�
syms t; % ������ű��� t

num=PRF*t_accumulate; % �ۻ�����
tau=t_accumulate*duty_ratio/num; % ÿ�������ʱ�䳤��
% ������ʱ��
tau_interval=t_accumulate/PRF; % ������ʱ��
%���β�
A=@(t) t<tau & t>=0;
for i=-num/2:num/2-1
    lower= i*tau;
    upper= (i+1)*tau;
    A(t)=A(t) | (t>i*tau_interval & t<i*tau_interval+tau);
end


% ���Ե�Ƶ���崮�źţ�����λ��
x=@(t) A(t).*exp(1i*(pi*k*t.^2+2*pi*f_c*t));

x_sample=x(t_sample); % �����ź�

