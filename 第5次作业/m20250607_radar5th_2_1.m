
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
R_list=[0,1e3,1.5e3, 2e3,2.5e3,3e3, 4e3, 8e3]; % Ŀ������б�
%R2_target=2e5; % Ŀ����루�ף�
c= 3e8; % ���٣���/�룩
tau= 2 * R_target / c; % �ز��źŵ���ʱ���룩

% ����
n = @(t) 0.1 * randn(size(t)); % ��һ��Ŀ�������

% ����ز��ź�
s=@(t) x(t)+n(t);

% ͨ��ƥ���˲���
y=conv(s(t_sample),fliplr(conj(x(t_sample))));

y_max= max(abs(y)); % ��һ��Ŀ������ֵ
% 3dB ��
y_3db = 20*log10(0.707); % 3dB��
% ���ƽ��
f1=figure(1);

title('ƥ���˲������');
xlabel('ʱ�� (��)');
ylabel('����');
%����3dB��
plot(t_sample+tau, y_3db*ones(size(t_sample)), 'k--', 'DisplayName', '3dB��');
hold on;


for index=1:1:length(R_list)
    R1_target =R_list(index);
    t1= 2 * R1_target / c;%���ʱ���
    s1 =@(t) x(t+t1)+n(t+t1);
    y1=conv(s1(t_sample),fliplr(conj(x(t_sample))));
    y1_db=20*log10(abs(y1)./y_max); % Ŀ���dBֵ
    t1_sample= t_new  + t1 + tau; % Ŀ���ʱ������
    plot(t1_sample,y1_db,'DisplayName', ['����' ,num2str(R1_target)]);
    hold on;
end
legend('show');
hold off;