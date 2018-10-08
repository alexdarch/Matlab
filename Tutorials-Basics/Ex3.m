clear all
close all

%question 3 exponential distribution
figure(1)
Num = 10000;
x=rand(Num,1);         %set up uniform distribution
w=0:0.01:7;
T = -log(1-x);

subplot(211),
hold on
U = histogram(T, 20, 'Normalization', 'pdf');  %plot histogram estimate of pdf
plot(w, exp(-w))                               %overlay exponential distribution
hold off

title('Histogram Estimate of CDF Method on an Exponential distribution'); %format plots
xlabel('y');
ylabel('Normalised Probability (p.d.f.)');
axis([0 7 0 inf]);
txt = sprintf('Histogram Estimate with %i Samples', Num);
legend(txt, 'Exponential pdf, f_Y(y)')


subplot(212),
width = 0.01;
hold on
ksdensity(T, 'width', width)    %plot smoothed kernel
plot(w, exp(-w))                %overlay exponential distribution
hold off

title('Histogram Estimate of CDF Method on an Exponential distribution'); %format plots
xlabel('y');
ylabel('Normalised Probability (p.d.f.)');
axis([-inf 7 0 inf]);
txt = sprintf('Kernel Density Estimate with width: %.2f', width);
legend(txt, 'Exponential pdf, f_Y(y)')



%FTR plot of mean squared error
N=200; spacing = 1;
ux=rand(N, 1);              %set up uniform distribution with N random values
Random_exp = -log(1-ux);            %get random samples

n = 1:1:N;
one_vector = ones(1, N);    %allows calculation of cumulative 1/n vector

MC_mean = (one_vector./n)'.*cumsum(Random_exp);     %monte carlo mean u(n) = (1/n)*sum(random variables up to n)
                                                    %note ' = transpose, and this is used because otherwise dimensions arent equal
%var_exp = (one_vector./n)'.*cumsum(Random_exp.^2) - MC_mean.^2; 
error = (MC_mean' - 1).^2;                          %the arithmetic (monte carlo) mean minus the theoretical mean (1)

figure(2)
hold on
plot(n, one_vector./n, 'b')
plot(n, error, 'r')
hold off

title('Plot of the Squared Error as N increases for the Exponential Distribution with \mu = 1, \sigma^2 = 1'); %format plots
xlabel('N');
ylabel('Mean Squared Error');
legend('Expected Error = \sigma^2/N', 'Cumulative Mean Squared Error')
axis([0 N 0 0.5]);




