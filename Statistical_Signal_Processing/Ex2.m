close all
clear all

%question 2a y=ax+b
Num_a = 10000; a = 0.5; b = 1;
E_xa = linspace(-10, 10, 1000);  %x components of expected transformed pdf
f_xa=randn(Num_a,1);            %generates Num_a normally distributed random numbers in a 1xNum_a column vector

T_ya= a*f_xa + b;                %transforms random numbers in 'x' to 'y'
E_fy = (1/sqrt(2*pi))*exp(-0.5*((E_xa-b)/a).^2)/a;   %expected transformed distribution

figure(1)   %plot graphs
hold on
H_a= histogram(T_ya, 20, 'normalization', 'pdf');    %normalised histogram of tranformed distribution
plot(E_xa, E_fy, 'g')   %plot transformed distribution
plot(E_xa, (1/sqrt(2*pi))*exp(-0.5*((E_xa)).^2), 'b') %plot distribution before transformation
hold off

title('Gaussian Distribution Transformed by y = ax+b'); %format plots
xlabel('x');
ylabel('Normalised Probability (p.d.f.)');
axis([-4 4 0 inf]);
txt = sprintf('Histogram Estimate with %i Samples', Num_a);
legend(txt, 'Transformed pdf, f_Y(y) (letting y=x)', 'Pre-transformation pdf, f_X(x)')

clear all





%question 2b y=x^2
Num_b = 10000;
E_xb = linspace(-10, 10, Num_b); %x components of expected transformed pdf
f_xb=randn(Num_b,1);    %generates Num_b normally distributed random numbers in a 1xNum_b column vector

T_yb = f_xb.^2;                                 %transforms random numbers in 'x' to 'y'
E_fy = (1./sqrt(2*pi*E_xb)).*exp(-0.5*E_xb);    %expected transformed distribution

figure(2)
hold on
H_b = histogram(T_yb, 20, 'Normalization', 'pdf', 'BinLimits', [-1 4]);      %plot normalised histogram plot
plot(E_xb, E_fy, 'g')                                   %plot transformed distribution
plot(E_xb, (1/sqrt(2*pi))*exp(-0.5*((E_xb)).^2), 'b')   %plot distribution before transformation
hold off

title('Gaussian Distribution Transformed by y = x^2'); %format plots
xlabel('x');
ylabel('Normalised Probability (pdf)');
axis([-2 5 0 2]);
txt = sprintf('Histogram Estimate with %i Samples', Num_b);
legend(txt, 'Transformed pdf, f_Y(y) (letting y = x)', 'Pre-transformation pdf, f_X(x)')

clear all




%FTR (c) y = sin(x)
Num_c = 10000; Uni = 1/(2*pi).*ones(Num_c, 1);

E_xc = linspace(0, 1/2*pi, Num_c);%xcomponents of pre-transformed pdf
E_yc = linspace(-1, 1, Num_c);  %x components of expected transformed pdf

f_xc = (2*pi).*rand(Num_c,1);    %generates Num_c uniformly distributed random numbers between 0&2*pi
                                   %so cdf=1 for interval 0<x<1/2pi (E_xc)
                                   
T_yc = sin(f_xc);                  %transforms random numbers in 'x' to 'y'
E_fy = (1/pi)./abs(sqrt(1-E_yc.^2));       %expected transformed distribution

figure(3)
hold on
H_c = histogram(T_yc, 40, 'Normalization', 'pdf', 'BinLimits', [-1 1]); %plots histogram estimation
plot(E_yc, E_fy, 'g')       %plots transformed distribution
plot(E_xc, Uni, 'b')        %plots original transformation
hold off

axis([-1 1 0 2.2]);
title('Uniform Distribution Transformed by y = sin(x)'); %format plots
xlabel('x');
ylabel('Normalised Probability (pdf)');
txt = sprintf('Histogram Estimate with %i Samples', Num_c);
legend(txt, 'Transformed pdf, f_Y(y) (letting y=x)', 'Pre-transformation pdf, f_X(x)')

clear all





%FTR (d) y = min(sin(x), 0.7)
Num_d = 10000; Uni = 1/(2*pi).*ones(Num_d, 1);

E_xd = linspace(0, 1/2*pi, Num_d);%xcomponents of pre-transformed pdf
E_yd = linspace(-1, 0.7, Num_d);  %x components of expected transformed pdf

f_xd = (2*pi).*rand(Num_d, 1);    %generates Num_c uniformly distributed random numbers between 0&2*pi           
                                   
T_yd = min(sin(f_xd), 0.7);                  %transforms random numbers in 'x' to 'y'
E_fy = (1/pi)./abs(sqrt(1-E_yd.^2));       %expected transformed distribution

figure(4)
hold on
H_c = histogram(T_yd, 40, 'Normalization', 'pdf', 'BinLimits', [-1 1]); %plots histogram estimation
plot(E_yd, E_fy, 'g')       %plots transformed distribution
hold off

axis([-1 1 0 6]);
title('Uniform Distribution Transformed by y = sin(x)'); %format plots
xlabel('x');
ylabel('Normalised Probability (pdf)');
txt = sprintf('Histogram Estimate with %i Samples', Num_d);
legend(txt, 'Transformed pdf, f_Y(y) (letting y=x)')
