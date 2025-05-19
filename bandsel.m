%% ============================= bandsel =============================
function [x, y] = bandsel(theta0, tau, n1)
% [x, y] = bandsel([57.98, 46.43]', 0.5, 10)

format short g
close ALL
rng(300); % For reproducibility

H = [0.25 1 3 10 20 40 50 80 100 200 300 400];
nh = length(H);

initheta = theta0  + 3.*randn(2,1); % Initial theta to run fminsearch to compute population quantile
XL = 5; % Lower limit of design space
XU = 400; % Upper limit of design space
N = 100; % Total number of points in design space
X = linspace(XL, XU, N); 

sigma1 = ones(1,N);
sigma2 = (X + X(N)).^2;
sigma3 = X + X(N);
sigma4 = (X + (0.25*X(N))).^2;

% Choose the true sigma
sigma = sigma2; 

% Initial design: Equispaced design
ix = linspace(XL, XU, n1)';
for j = 1:n1
    ix(j) = closest(X', ix(j));
end

%Data generation
m = 20; 
F0all = (theta0(1)*X')./(theta0(2)+X'); % Generated using theta0
%%Derivative of F w.r.t theta

Xv = zeros(N,m);
Yv = zeros(N,m); %Generate repeated responses for population

for j = 1:m
 Xv(:,j) = X';
end 

% Standardized Variance function
Tsigma = sqrt(sigma'/mean(sigma));

% Generating errors from standard normal distribution
a = 2; %sigma epsilon
b = norminv(tau);

for j = 1:m
 errors = (a*Tsigma).*randn(N,1) - Tsigma*a*b;
 Yv(:,j) = F0all + errors;
end


eIND = (1:N)';
IX = zeros(N, 2);
IX(:,1) = eIND;
IX(:,2) = X';

Ind = 1:N;

x = ix;
n1 = length(x);
n = n1;
sInd = Ind(ismember(X',x)); % Initial sample indices
y = zeros(n,1);
for j = 1:n
    y(j) = datasample(Yv(sInd(j),:), 1);
end

% Quntile objective function that is required to minimize
fminfun = @(theta) sum(((y - theta(1)*x./(theta(2) + x))).*(tau - vindicator(y - theta(1)*x./(theta(2) + x))));
thetahat = fminsearch(fminfun, initheta); % estimates of theta

Sh = zeros(n, nh);

for i = 1:nh
[regdata] = bandscales(x, y, N, IX, H(i), thetahat);
Sh(:,i) = regdata(:,2);
end

for i = 1:nh
   subplot(4, 3, i);
   scatter(x, Sh(:,i));
   hold on  
   labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)'];
   xlabel(labels((3*(i-1)+1):(3*(i-1)+3)));
   ylabel('S_n(x)');
end

%% ============================= bandscales =============================
function [regdata] = bandscales(x, y, N, IX, h, thetahat)

sdx = zeros(N, 1);
np = zeros(N, 1);

for j = 1:N
np(j) = sum(ismember(x, IX(j, 2)));
end

nonz = sum(np > 0);
regdata = zeros(nonz, 2);
k = 0;
for j = 1:N
    if np(j) > 0
        k = k + 1;
%Estimation of scale functions
KW = Scale(IX(j, 2), x, h);

merror = y - Fun(x, thetahat); % Model errors
sdx(j) = sqrt(sum(KW.*(merror.^2))); %Estimates of the scale function
regdata(k ,1) = IX(j, 2); %Regression x variable
regdata(k ,2) = sdx(j);
    end
end

%% ============================= closest =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xmin = closest(X,y)
% finds the closest member of X to y.
dist = abs(X-y);
d0 = min(dist);
xmin = X(dist==d0);
if length(xmin) > 1
    xmin = xmin(1);
end

%% ============================= Fun =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = Fun(x, theta)  
    % Forms the vector of nonlinear, Michaelis-Menten,  responses
    % theta is parameter (2 x 1)
    % x is a design vector (n x 1), where n is the sample size
    % out (n x 1) 
out = theta(1)*x./(theta(2) + x);

%% ============================= indicator =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outv = indicator(r) 
% Indicator function
if r < 0
   outv = 1;
else
   outv = 0;
end

%% ============================= Scale =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = Scale(xv, x, h) 
    % xv is the input (1 x 1)
    % x is the sample design points (n x 1)
    % y is the sample responses (n x 1)
    % thetahat is the tau^th quantile estimates using sample (p x 1)
    % h is the bandwidth 
    kweight = exp((-0.5)*(((xv - x)./h).^2)); % Standard Gaussian kernel
    result = kweight./sum(kweight); % Standardized weights
 
    %% ============================= vindicator =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outv = vindicator(vr)
% vector indicator function 
% vr is vector (n x 1)
n = length(vr);
outv = zeros(n,1);
for i = 1:n
   if vr(i) < 0
   outv(i) = 1;
   else
   outv(i) = 0;
   end 
end
