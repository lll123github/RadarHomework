sample_num=2400;
theta=-pi:2*pi/sample_num:pi;
x_axis=-180:360/sample_num:180;
lambda=0.1;
N=20;
% 固定theta_0=0，比较不同的d带来的影响
f1=figure(1);
subplot_index=0;
for magnification=[0.8 1 1.2]
    d=magnification*lambda;
    subplot_index=subplot_index+1;
    tmp_phi= 2*pi*d/lambda*sin(theta);
    E_a_theta=zeros(1,length(tmp_phi));
    for index=1:length(tmp_phi)
        if abs(tmp_phi(index))<10e-5 || abs(tmp_phi(index)-2*pi)<10e-5
            E_a_theta(index)=N;
        else
            E_a_theta(index)=abs(sin(N*tmp_phi(index)/2)./sin(tmp_phi(index)/2));
        end
    end
    G_a_theta=E_a_theta.^2./(N^2);
    G_a_theta_db=10*log10(G_a_theta);
    figure(f1);
    subplot(3,1,subplot_index);
    plot(x_axis, G_a_theta_db);
    ylim([-70 0]);
    xlim([-180 180]);
    xlabel('\theta(度)');
    ylabel('G_a(\theta) (dB)');
    title(['d=', num2str(magnification), '\lambda']);
    hold on;
end

% 画出不同观测角度的方向图 %固定d
f2=figure(2);
subplot_index=0;
for theta_0=[0 pi/12 pi/6 pi/4 pi/3]
    
    d=0.8*lambda;
    tmp_phi= 2*pi*d/lambda*(sin(theta)-sin(theta_0));
    E_a_theta=zeros(1,length(tmp_phi));
    for index=1:length(tmp_phi)
        if abs(tmp_phi(index))<10e-5 || abs(tmp_phi(index)-2*pi)<10e-5
            E_a_theta(index)=N;
        else
            E_a_theta(index)=abs(sin(N*tmp_phi(index)/2)./sin(tmp_phi(index)/2));
        end
    end
    G_a_theta=E_a_theta.^2./(N^2);
    G_a_theta_db=10*log10(G_a_theta);
    if theta_0==0
        G_a_theta_db_origin=G_a_theta_db;
    end
    figure(f2);
    subplot_index=subplot_index+1;
    subplot(5,1,subplot_index);
    plot(x_axis+theta_0*180/pi, G_a_theta_db_origin);
    hold on;
    plot(x_axis, G_a_theta_db);
    ylim([-70 0]);
    xlim([-180 180]);
    xlabel('\theta(度)');
    ylabel('G_a(\theta) (dB)');
    title(['\theta_0=', num2str(theta_0*360/2/pi), '度']);
    hold on;

    % figure(f3);
    % subplot(4,1,subplot_index);
    % plot(x_axis, circshift(G_a_theta_db_origin, theta_0/2/pi*sample_num));
    % hold on;
    % plot(x_axis, G_a_theta_db);
    % ylim([-70 0]);
    % xlim([-180 180]);
    % xlabel('\theta(度)');
    % ylabel('G_a(\theta) (dB)');
    % title(['\theta_0=', num2str(theta_0*360/2/pi), '度']);
    % hold on;
end

%固定phi
f3=figure(3);
subplot_index=0;
for theta_0=[0 pi/12 pi/6 pi/4 pi/3]
    subplot_index=subplot_index+1;
    d=0.5*lambda;
    phi= 2*pi*d/lambda*sin(theta_0);
    tmp_phi= 2*pi*d/lambda*sin(theta)-phi;
    E_a_theta=zeros(1,length(tmp_phi));
    for index=1:length(tmp_phi)
        if abs(tmp_phi(index))<10e-5 || abs(tmp_phi(index)-2*pi)<10e-5
            E_a_theta(index)=N;
        else
            E_a_theta(index)=abs(sin(N*tmp_phi(index)/2)./sin(tmp_phi(index)/2));
        end
    end
    G_a_theta=E_a_theta.^2./(N^2);
    G_a_theta_db=10*log10(G_a_theta);

    lambda1=0.95*lambda;
    tmp_phi1= 2*pi*d/lambda1*sin(theta)-phi;
    E_a_theta1=zeros(1,length(tmp_phi1));
    for index=1:length(tmp_phi1)
        if abs(tmp_phi1(index))<10e-5 || abs(tmp_phi1(index)-2*pi)<10e-5
            E_a_theta1(index)=N;
        else
            E_a_theta1(index)=abs(sin(N*tmp_phi1(index)/2)./sin(tmp_phi1(index)/2));
        end
    end
    G_a_theta1=E_a_theta1.^2./(N^2);
    G_a_theta_db1=10*log10(G_a_theta1);
    figure(f3);
    subplot(4,1,subplot_index);
    plot(x_axis, G_a_theta_db);
    hold on;
    plot(x_axis, G_a_theta_db1);
    ylim([-70 0]);
    xlim([-180 180]);
    xlabel('\theta(度)');
    ylabel('G_a(\theta) (dB)');
    title(['\theta_0=', num2str(theta_0*360/2/pi), '度']);
    hold on;
end

