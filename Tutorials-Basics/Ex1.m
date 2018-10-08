clear all
close all
%% 
%smoothed uniform distribution:
Num_U = 1000;                       %number random numbers generated
U1=rand(Num_U,1);                   %generate Num_U uniformily distributed random numbers in a 1xNum_U column vector
E_Ux = linspace(0, 1, Num_U);       %x components of expected uniform distribution
E_Uy = ones(Num_U, 1);              %expected distribution

[K_U, K_x] = ksdensity(U1, 'width', 0.1); %returns smoothed distribution K_U = fn(K_x)
                                          %sets random distribution to smooth as U1, with kernel width = 0.1

%plot histogramed uniform distribution
figure(1)
j = 1;
for Num_Uh = [100 1000 10000]
subplot(3, 1, j),                 %generate three vertically aligned subplots
U2 = rand(Num_Uh, 1);             %generate Num_Uh uniformly distributed random numbers

hold on
H_U = histogram(U2, 10);                       %stores bin information in a struct and plots 10 bins
plot(E_Ux, E_Uy.*(Num_Uh/H_U.NumBins) , 'g')   %plot expected distribution

i=1:1:10;                                    %calculate mean and std for each of the 10 bars
mean(i) = Num_Uh/H_U.NumBins;          
std(i) = sqrt(Num_Uh*H_U.BinWidth*(1-H_U.BinWidth)); %note division by average numel per bin
errorbar((i/10)-0.05, mean, 3*std)           %plot 3*std errorbars for all 10 bars

axis([0 1 0 inf]);
ttl = sprintf('Number of Samples: %i', Num_Uh);
title(ttl);
ylabel('Histogram Count');
xlabel('x');
hold off

j=j+1;  %next subplot
end

%comparison of kernel density and histogram plot
figure(2)
subplot(2, 1, 1),
title('Kernel Density Method');
ylabel('f_X(x)');
xlabel('x');
axis([0 1 0 2]);

hold on
plot(K_x, K_U, 'b')
plot(E_Ux, E_Uy, 'g')
hold off


subplot(2, 1, 2),
title('Histogram Method');
ylabel('f_X(x)');
xlabel('x');
axis([0 1 0 2]);

hold on
H_U = histogram(U1, 10, 'Normalization', 'pdf');
plot(E_Ux, E_Uy, 'g')
hold off
%% 



clear all
%smoothed gaussian distribution
Num_G = 10000;                             %number random numbers generated
G1 = randn(Num_G,1);                       %generate Num_U normally distributed random numbers in a 1xNum_G column vector
E_Gx = linspace(-5, 5, Num_G);             %x components of expected gaussian distribution
E_Gy = (1/sqrt(2*pi)).*exp(-0.5*(E_Gx.^2)); %y components of expected gaussian distribution

[K_G, K_x] = ksdensity(G1, 'width', 0.1);  %returns smoothed distribution K_G = f(K_x)
                                           %sets random distribution to smooth as G1, with kernel width = 0.1

%plot histogrammed normal distribution:
figure(3)
j = 1;
for Num_Gh = [100 1000 10000]
subplot(3, 1, j),           %generate three vertically aligned subplots
G2 = randn(Num_Gh, 1);      %generate Num_Gh uniformly distributed random numbers

hold on
H_G = histogram(G2, 20, 'Normalization', 'count');  %gives a histogram where the heights are the number of elements per bin
plot(E_Gx, E_Gy*H_G.BinWidth*Num_Gh, 'g')           %plot scaled gaussian distribution (use bin width and not num bins because delta =! 1/J anymore)222
start = H_G.BinLimits(1);                   %define the x coordinate of the edge of the far left bin

mean = zeros(1, 20);  std = zeros(1, 20); prob = zeros(1, 20);

    for i=1:1:20 %calculate mean and std for each of the 20 bars
        binStart = start+((i-1)*H_G.BinWidth);  %left edge coordinate of bin i
        binEnd = start+(i*H_G.BinWidth);        %right edge coordinate of bin i
        q = normcdf([binStart binEnd]);         %returns the cdf at x=binStart and x=binEnd of a gaussian
        prob(i) = q(2) - q(1);                  %returns p(binStart < x < binEnd)

       
        mean(i) = prob(i)*Num_Gh;               %for un-normalised histogram mean, multiply by Num_Gh
                                                %note, this is the theoretical mean, not the experimental one
        std(i) = sqrt(Num_Gh*prob(i)*(1-prob(i)));
        errorbar(binStart+(H_G.BinWidth/2), mean(i), 3*std(i))          %plot 3*std errorbars for all 20 bars
        
    end
ttl = sprintf('Number of Samples: %i', Num_Gh);
title(ttl);
ylabel('Histogram Count');
xlabel('x');
axis([-4 4 0 inf])      %sets ymin=0 but allows for auto ymax

hold off
j=j+1;
end

%comparison of kernel density and histogram plot
figure(4)
subplot(2, 1, 1),
title('Kernel Density Method');
ylabel('f_X(x)');
xlabel('x');
axis([-4 4 0 inf]);

hold on                                        
plot(K_x, K_G, 'b')
plot(E_Gx, E_Gy, 'g')
hold off

subplot(2, 1, 2),
title('Histogram Method');
ylabel('f_X(x)');
xlabel('x');
axis([-4 4 0 inf]);

hold on
H_G = histogram(G1, 20, 'Normalization', 'pdf');
plot(E_Gx, E_Gy, 'g')
hold off

clear all
figure(5) %plot of bin probabilities vs variance

p = linspace(0,1, 200);
one = ones(1, 200);
var = p.*(one - p);
p_max = 1/(2*pi);

hold on
plot(p, var, 'b')
plot([0 p_max p_max], [p_max*(1-p_max) p_max*(1-p_max) 0], 'g.-')
hold off

xlabel('p_j')
ylabel('\sigma_j^2 = p_j(1-p_j)')
title('plot of \sigma_j^2 vs p_j')





