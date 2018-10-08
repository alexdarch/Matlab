close all
clear all

%question 4

%plot differing values of alpha and beta to compare effects
figure(1)
j=1;
for alpha=[0.5 1.5]     %Plot different values of alpha
    for beta=-1:0.5:1   %Plot different values of beta

    Num = 100000;

    b = atan(beta*tan(pi*alpha/2))/alpha;          %difficult density algebra
    s =(1+(beta*tan(pi*alpha/2))^2)^(1/(2*alpha));
    
    U_rand = -pi/2 + pi*rand(Num, 1);   %generate random numbers between -pi/2 < x < pi/2
    E_rand = exprnd(1, Num, 1);         %generate exponential distribution with mean = 1

    X = s*(sin(alpha*(U_rand+b))./(cos(U_rand)).^(1/alpha)).*((cos(U_rand - alpha*(U_rand+b)))./E_rand).^((1-alpha)/alpha);     %generate 'difficult' density
    
    X(abs(X) > 30) = 0;             %replace values of X > 30 with zero (removes outliers)
    graph = sprintf('Alpha = %0.1f, Beta = %0.1f', alpha, beta); 
    
    subplot(2, 5, j),               %plot histograms
    histogram(X, 'BinWidth', 1, 'normalization', 'pdf')
    
    title(graph)
    xlabel('x')
    ylabel('Difficult pdf, f_X(x)')
    axis([-20 20 0 0.5])       %graph size limits so comparable
    
    j = j+1;
  
    end
end

clear all





%relationship between alpha and the variance
figure(2)
Num = 100000;   beta = 0;

i=1;
for alpha=0.01:0.01:2     %Plot different values of alpha

%same calculation as fisrt part
b = atan(beta*tan(pi*alpha/2))/alpha;          %difficult density algebra
s =(1+(beta*tan(pi*alpha/2))^2)^(1/(2*alpha));

U_rand = -pi/2 + pi*rand(Num, 1);   %generate random numbers between -pi/2 < x < pi/2
E_rand = exprnd(1, Num, 1);         %generate exponential distribution with mean = 1

X = s*(sin(alpha*(U_rand+b))./(cos(U_rand)).^(1/alpha)).*((cos(U_rand - alpha*(U_rand+b)))./E_rand).^((1-alpha)/alpha);     %generate 'difficult' density
X_abs = abs(X);
X(abs(X) > 30) = 0;             %replace values of X > 30 with zero (removes outliers) 
    
  

%Monte-Carlo variance estimator    
q = 1:1:Num;
one_vector = ones(1, Num);    %allows calculation of cumulative 1/n vector
MC_mean = (one_vector./q)'.*cumsum(X);     %monte carlo mean u(n) = (1/n)*sum(random variables up to n)
MC_var = (one_vector./q)'.*cumsum(X.^2) - MC_mean.^2;   %monte carlo variance estimate

variance_est(i) = MC_var(numel(MC_var)-1);
    
i = i+1;    %next value
end

alpha=0.01:0.01:2;
hold on
plot(alpha, variance_est, 'r')
plot(alpha, 238.8.*alpha.*exp(-2.5437.*alpha), 'b')
hold off


title('Plot of Variance vs Stability Parameter')
legend('Monte-Carlo Variance Estimate', 'Plot of 238.8\alphaexp(-2.5437\alpha)')
xlabel('\alpha')
ylabel('Monte-Carlo Variance estimate')
axis([0 inf 0 inf]) 




clear all








%plot tail probabilities
figure(3)
Num = 100000;   beta = 0; alpha = [0.5 1 1.5 2];


for i=1:1:numel(alpha)    %Plot different values of alpha

%same calculation as last part
b = atan(beta*tan(pi*alpha(i)/2))/alpha(i);          %difficult density algebra
s =(1+(beta*tan(pi*alpha(i)/2))^2)^(1/(2*alpha(i)));

U_rand = -pi/2 + pi*rand(Num, 1);   %generate random numbers between -pi/2 < x < pi/2
E_rand = exprnd(1, Num, 1);         %generate exponential distribution with mean = 1

X_abs = abs(s*(sin(alpha(i)*(U_rand+b))./(cos(U_rand)).^(1/alpha(i))).*((cos(U_rand - alpha(i)*(U_rand+b)))./E_rand).^((1-alpha(i))/alpha(i)));     %generate 'difficult' density
X_abs(X_abs > 30) = 0;             %replace values of X > 30 with zero (removes outliers) 
    

    k=1;    %calculate tail probabilities
    for t = 0:0.1:10*sqrt(2);
        m=1;
        n = 0;
        while(m < Num+1)    %limit size of tail probability
            if (X_abs(m)>t)     %check if the mth element is greater than the tail probability
                n = n+1;        %if it is then increase the number of elements greater than t
            end
            m = m+1;        %check next element
        end
        tailP(k) = n/Num;   %(num_random_values > t) normalised = tail probability
        gausTailP(k) = (1-normcdf(t))*2;    %gaussian tail prob for comparison
        k=k+1;              %next tail prob for next value of t
    end
    
subplot(numel(alpha), 1, i),   
t = 0:0.1:10*sqrt(2);
hold on
plot(t./sqrt(2), tailP, 'r')
plot(t, gausTailP, 'b')
hold off   


graph = sprintf('Alpha = %0.1f, Beta = %0.1f', alpha(i), beta); 
title(graph)
xlabel('t')
ylabel('p(|X|>t))')
axis([0 10 0 1])       %graph size limits so comparable

tail_alpha(:, i) = tailP;   %store tailP in columns of differing alpha for use in next section.

end



%plot tail probabilities vs gamma estimates
figure(4)

for r = 1:1:numel(alpha)
    
subplot(numel(alpha), 1, r),    
hold on
plot(t./sqrt(2), tail_alpha(:, r), 'r')
plot(t, 0.35.*(t.^(-alpha(r))), 'b')
hold off   


graph = sprintf('Alpha = %0.1f, Beta = %0.1f', alpha(r), beta); 
title(graph)
xlabel('t')
ylabel('p(|X|>t))')
axis([0 10 0 1])       %graph size limits so comparable
    
end
    
