tau=2;%���Ե�Ƶ�ź�ʱ�䳤��
t_all=4;%�ź���ʱ��
k=50;%б��
B=k*tau;%����
% f_c=300;%����Ƶ��
fs=50*B;%������
ts=1/fs;%�������
N=floor(t_all*fs)+1;%��������(�����˵�)

% �������ʱ���
t= [-(N-1)/2:(N-1)/2]*ts;
% �� -1 �� 1������ 2001 ����
% ���֮���ʱ���
t_new=[-(N-1):N-1] * ts ;
%Ƶ��Ļ��Ʒ�Χ��Ƶ�ʷֱ���Ӧ����fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % Ƶ�ʷ�Χ�� -N/2 �� N/2

%���β�
A=t<tau/2&t>-tau/2;
% ���Ե�Ƶ�źţ�����λ��
x=A.*exp(1i*(pi*k*t.^2));

I=A.*cos(pi*k*t.^2);
Q=A.*sin(pi*k*t.^2);
S=I+1i.*Q;
match_filter_S=fliplr(conj(S));
y=conv(S,match_filter_S);

center_index = floor(length(y)/2);
range = 100;
t_index_start = center_index - range;
t_index_end = center_index + range;

f1=figure(1);
%����һ������µ�ƥ���˲����
plot(t_new(t_index_start:t_index_end),abs(y(t_index_start:t_index_end))); % ����ʵ������Χ�� -1 �� 1

f2=figure(2);
f3=figure(3);
f4=figure(4);
f5=figure(5);
f6=figure(6);
f7=figure(7);
f8=figure(8);
f9=figure(9);
f10=figure(10);

% ���Ȳ�ƽ��
for epsilon=[0,5,10]
    I_1=A.*(1+epsilon).*cos(pi*k*t.^2);
    Q_1=A.*sin(pi*k*t.^2);
    S_1=I_1+1i.*Q_1;
    y_1=conv(S_1,match_filter_S);
    figure(f2);
    plot(t_new(t_index_start:t_index_end),abs(y_1(t_index_start:t_index_end))); % ����ʵ������Χ�� -1 �� 1
    hold on;
    figure(f5);
    y_1_db=20*log10(abs(y_1)/max(abs(y_1))); % ת��ΪdB
    plot(t_new(t_index_start:t_index_end),y_1_db(t_index_start:t_index_end)); % ����ʵ������Χ�� -1 �� 1
    hold on;
    figure(f8);
    %������λ���
    y_1_phase=angle(y_1);
    plot(t_new(t_index_start:t_index_end),y_1_phase(t_index_start:t_index_end)); % ������λ
    hold on;

end
figure(f2);
xlabel('t');
ylabel('|y(t)|');
title('���Ȳ�ƽ�⣨��ƽ��');
legend('���Ȳ�Ϊ0','���Ȳ�Ϊ5','���Ȳ�Ϊ10');
figure(f5);
xlabel('t');
ylabel('20*lg(abs(y(t)))');
title('���Ȳ�ƽ�⣨dB��');
legend('���Ȳ�Ϊ0','���Ȳ�Ϊ5','���Ȳ�Ϊ10');
figure(f8);
xlabel('t');
ylabel('phase(y(t))');
title('���Ȳ�ƽ�⣨��λ��');
legend('���Ȳ�Ϊ0','���Ȳ�Ϊ5','���Ȳ�Ϊ10');


%��λ������
for delta_theta=[0,pi/4,pi/2]
    I_2=A.*cos(pi*k*t.^2-delta_theta);
    Q_2=A.*sin(pi*k*t.^2);
    S_2=I_2+1i.*Q_2;
    y_2=conv(S_2,match_filter_S);
    figure(f3);
    plot(t_new(t_index_start:t_index_end),abs(y_2(t_index_start:t_index_end))); 
    % ����ʵ������Χ�� -1 �� 1
    hold on;
    figure(f6);
    y_2_db=20*log10(abs(y_2)/max(abs(y_2))); % ת��ΪdB
    plot(t_new(t_index_start:t_index_end),y_2_db(t_index_start:t_index_end)); % ����ʵ������Χ�� -1 �� 1
    hold on;
    figure(f9);
    %������λ���
    y_2_phase=angle(y_2);
    plot(t_new(t_index_start:t_index_end),y_2_phase(t_index_start:t_index_end)); % ������λ
    hold on;
end
figure(f3);
xlabel('t');
ylabel('|y(t)|');
title('��λ������');
legend('��λ���Ϊ0','��λ���Ϊ\pi/2','��λ���Ϊ\pi');
figure(f6);
xlabel('t');
ylabel('20*lg(abs(y(t)))');
title('��λ��������dB��');
legend('��λ���Ϊ0','��λ���Ϊ\pi/4','��λ���Ϊ\pi/2');
figure(f9);
xlabel('t');
ylabel('phase(y(t))');
title('��λ����������λ��');
legend('��λ���Ϊ0','��λ���Ϊ\pi/4','��λ���Ϊ\pi/2');


% ֱ��ƫ��

for delta_k=[0,5,10]
    I_3=A.*(1).*cos(pi*k*t.^2)+delta_k;
    Q_3=A.*sin(pi*k*t.^2);
    S_3=I_3+1i.*Q_3;
    y_3=conv(S_3,match_filter_S);
    figure(f4);
    plot(t_new(t_index_start:t_index_end),abs(y_3(t_index_start:t_index_end))); 
    hold on;
    figure(f7);
    y_3_db=20*log10(abs(y_3)/max(abs(y_3))); % ת��ΪdB
    plot(t_new(t_index_start:t_index_end),y_3_db(t_index_start:t_index_end)); % ����ʵ������Χ�� -1 �� 1
    hold on;
    figure(f10);
    %������λ���
    y_3_phase=angle(y_3);
    plot(t_new(t_index_start:t_index_end),y_3_phase(t_index_start:t_index_end)); % ������λ
    hold on;
end
figure(f4);
xlabel('t');
ylabel('|y(t)|');
title('ֱ��ƫ��');
legend('ͬ�����ֱ��ƫ��Ϊ0','ͬ�����ֱ��ƫ��Ϊ5','ͬ�����ֱ��ƫ��Ϊ10');

figure(f7);
xlabel('t');
ylabel('20*lg(abs(y(t)))');
title('ֱ��ƫ�ã�dB��');
legend('ͬ�����ֱ��ƫ��Ϊ0','ͬ�����ֱ��ƫ��Ϊ5','ͬ�����ֱ��ƫ��Ϊ10');

figure(f10);
xlabel('t');
ylabel('phase(y(t))');
title('ֱ��ƫ�ã���λ��');
legend('ͬ�����ֱ��ƫ��Ϊ0','ͬ�����ֱ��ƫ��Ϊ5','ͬ�����ֱ��ƫ��Ϊ10');







