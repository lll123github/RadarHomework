
%�������Ե�Ƶ�ź�

tau=2e-5;%���Ե�Ƶ�ź�ʱ�䳤��
t_all=4e-5;%�ź���ʱ��
k=1e8;%б��
B=k*tau;%����
f_c=3e4;%����Ƶ��
fs=200*max(f_c,B);%������
ts=1/fs;%�������
N=floor(t_all*fs)+1;%��������(�����˵�)

% �������ʱ���
t_sample= [-(N-1)/2:(N-1)/2]*ts;
% �� -1 �� 1������ 2001 ����
% ���֮���ʱ���
t_new = [-(N-1):N-1] * ts ;
%Ƶ��Ļ��Ʒ�Χ��Ƶ�ʷֱ���Ӧ����fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % Ƶ�ʷ�Χ�� -N/2 �� N/2

syms t; % ������ű��� t
%���β�
A=@(t) t<tau/2&t>-tau/2;
% ���Ե�Ƶ�źţ�����λ��
x=@(t) A(t).*exp(1i*(pi*k*t.^2+2*pi*f_c*t));

x_sample=x(t_sample); % �����ź�

%�����ز��ź�

% Ŀ�����
R_target =1e5; % Ŀ����루�ף�
R2_target=R_target+4e3;
%R2_target=2e5; % Ŀ����루�ף�
c= 3e8; % ���٣���/�룩
t1= 2 * R_target / c; % �ز��źŵ���ʱ���룩
t2= 2 * R2_target / c; % �ز��źŵ���ʱ���룩

% ����
n1 = @(t) 0.1 * randn(size(t)); % ��һ��Ŀ�������
n2 = @(t) 0.1 * randn(size(t)); % �ڶ���Ŀ�������

% ����ز��ź�
s1=@(t) x(t)+n1(t);
s2=@(t) x(t+t2-t1)+n2(t+t2-t1);

% s1=@(t) x(t);
% s2=@(t) x(t+t2-t1);

s1_sample=s1(t_sample);
s2_sample=s2(t_sample);

% ͨ��ƥ���˲���
y1=conv(s1(t_sample),fliplr(conj(x(t_sample))));
y2=conv(s2(t_sample),fliplr(conj(x(t_sample))));

y1_max= max(abs(y1)); % ��һ��Ŀ������ֵ
y1_db= 20*log10(abs(y1)./y1_max); % ��һ��Ŀ���dBֵ
y2_db= 20*log10(abs(y2)./y1_max); % �ڶ���Ŀ���dBֵ

% �����Ӧʱ������
t1_sample= t_new  + t1; % ��һ��Ŀ���ʱ������
t2_sample= t_new  + t2; % �ڶ���Ŀ���ʱ������

% ���ƽ��
f1=figure(1);
% plot(t_sample,real(x(t_sample)));
% hold on;
% ����3dB ��
y_3db = 20*log10(0.707); % 3dB��

% ��һ����y1_db������ͬ��3dB�ߣ�x����t1_sampleһ��
plot(t1_sample, y_3db*ones(size(t1_sample)), 'k--', 'DisplayName', '3dB��');
hold on;
plot(t1_sample,abs(y1),'DisplayName', 'Ŀ��1');
hold on;
plot(t2_sample,abs(y2),'DisplayName', 'Ŀ��2');
title('ƥ���˲������');
xlabel('ʱ�� (��)');
ylabel('����');
% legend('Ŀ��1','Ŀ��2');
hold off;






