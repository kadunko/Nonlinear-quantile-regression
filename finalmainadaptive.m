function [AdaptiveTable] = finalmainadaptive(M, ds, plotwhich, theta0, eta, tau)
% M is the number of simulations
% ds is the desired sample size.
% plotwhich is the desired values of nu
% eta is the bound of model errors
% theta0 is true parameters that is required to generate responses in adaptive designs.
% tau is desired quantile level.
% [AdaptiveTable] = finalmainadaptive(1, 1000, [.2 .35 .65 .95], [57.98, 46.43]', 1, 0.5)

format short g
close ALL
rng(200); % For reproducibility

XL = 5; % Lower limit of design space
XU = 400; % Upper limit of design space
N = 100; % Total number of points in design space
X = linspace(XL, XU, N); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial design with 2 design points
%ix = [X(2), X(75)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial design: Equispaced design
n1 = 10;
ix = linspace(XL, XU, n1)';
for j = 1:n1
    ix(j) = closest(X', ix(j));
end

lpw = length(plotwhich);
nrows = 4; % the number of scale functions
MRCxpoints = zeros(ds, lpw, nrows, M);
finalxpoints = zeros(ds, lpw, nrows);
finaldesign = zeros(N, lpw, nrows);
FinalAtrueloss = zeros(lpw, nrows, M);
AverageAtrueloss = zeros(lpw, nrows);
finalRESIGMA = zeros(N,lpw, nrows, M);
AverageRESIGMA = zeros(N, lpw, nrows);

for m = 1:M
    [RCxpoints, Atrueloss, RTsigma, RESIGMA] = mainadaptive(ds, plotwhich, theta0, eta, tau, ix, XL, XU, N);
    MRCxpoints(:,:,:,m) = RCxpoints;
    FinalAtrueloss(:,:,m) = Atrueloss;
    FinalRESIGMA(:,:,:,m) = RESIGMA;
end

% Computing averaged true losses
for count = 1:lpw
    for row = 1:nrows
        AverageAtrueloss(count, row) = mean(FinalAtrueloss(count, row,:));
    end
end

% Adaptive Loss Table
rowNames = {'nu = 0.20', 'nu = 0.35', 'nu = 0.65', 'nu = 0.95'};
colNames = {'Sigma1', 'Sigma2', 'Sigma3', 'Sigma4'};
AdaptiveTable = array2table(AverageAtrueloss,'RowNames',rowNames,'VariableNames',colNames);

for i = 1:N
    for count = 1:lpw
        for row = 1:nrows
            AverageRESIGMA(i, count, row) = mean(FinalRESIGMA(i, count, row,:));
        end
    end
end


for s = 1:ds
    for count = 1:lpw
        for row = 1:nrows
            meanxpoint = mean(MRCxpoints(s,count,row,:));
            finalxpoints(s, count, row) = closest(X,meanxpoint);
        end
    end
end

for s = 1:ds
    for count = 1:lpw
        for row = 1:nrows
            for i = 1:N
                if isequal(finalxpoints(s, count, row), X(i))
                   finaldesign(i, count, row) = finaldesign(i, count, row) + 1;
                end
            end
        end
    end
end
    
finaldesign = finaldesign/ds;

figure(1) %% Plots of designs
for row = 1:nrows
    for count = 1:lpw
    subplot(nrows, lpw, lpw*(row - 1) + count);
    mm1 = 1.2*max(finaldesign(:, count, row));
    mm2 = 1.5*max(finaldesign(:, count, row));
    mm = min(mm1,mm2);
    bar(X, finaldesign(:, count, row), 'EdgeColor', 'b', 'FaceColor', 'none');
    axis([1 400 -.05*mm mm1]);
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)'];
    group = 4*(row - 1) + count;
    xlabel(labels((3*(group-1)+1):(3*(group-1)+3)));
	ylabel('\xi(x)'); 
    end
end

Maxt = 0;
for count = 1:lpw
       Maxc = max(RTsigma(:,1,count));
       Maxt = max(Maxc, Maxt);
end

Maxe = 0;
for row = 1:nrows
   for count = 1:lpw
       Maxcr = max(AverageRESIGMA(:, count, row));
       Maxe = max(Maxcr, Maxe);
   end
end

Maxv = max(Maxt, Maxe) + 1;

figure(3) %% Plots of designs
for row = 1:nrows
for count = 1:lpw
    subplot(nrows, lpw, lpw*(row - 1) + count);
    SSigma = RTsigma(:,:,row).^2;
    EstSigma = AverageRESIGMA(:,count, row).^2;
    
    plot(X', SSigma);
    hold on;
    plot(X', EstSigma);
    ylim([0 Maxv]);
    xlim([XL XU]);
    legend('T', 'E','Location','northwest');
    
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)'];
    
    group = 4*(row - 1) + count;
    xlabel(labels((3*(group-1)+1):(3*(group-1)+3)));
    ylabel('\sigma(x)^2');
    
end
end

%% ====================== Main function to choose design points ========================
function [RCxpoints, Atrueloss, RTsigma, RESIGMA] = mainadaptive(ds, plotwhich, theta0, eta, tau, ix, XL, XU, N)

initheta = theta0  + 3.*randn(2,1); % Initial theta to run fminsearch to compute population quantile
X = linspace(XL, XU, N); % Design space
lpw = length(plotwhich);
Nl = 10000; %Preallocate size to save losses

%%Derivative of F w.r.t theta
Tfall = [X'./(theta0(2)+X'), -(theta0(1)*X')./((theta0(2)+X').^2)]; 
p = size(Tfall,2);
[TQ, ~] = qr(Tfall);
TQ1 = TQ(:,1:p);

% True variance functions
Rvar1 = ones(1,N);
Rvar2 = (X + X(N)).^2;
Rvar3 = X + X(N);
Rvar4 = ((2*X) + (0.5*X(N))).^2;

% Standardized sigma functions for satisfying the condition (2ii)
Rsigma1 = sqrt(Rvar1/mean(Rvar1)); 
Rsigma2 = sqrt(Rvar2/mean(Rvar2));
Rsigma3 = sqrt(Rvar3/mean(Rvar3));
Rsigma4 = sqrt(Rvar4/mean(Rvar4));

Realsigma = [Rsigma1; Rsigma2; Rsigma3; Rsigma4];

% Variance functions for adaptive designs
sigma1 = ones(1,N);
sigma2 = (X + X(N)).^2;
sigma3 = X + X(N);
sigma4 = (X + (0.25*X(N))).^2;
SIGMA = [sigma1; sigma2; sigma3; sigma4];
nrows = size(SIGMA, 1);

Atrueloss = zeros(lpw, nrows);

ROUT0 = zeros(N, lpw, nrows); % X versus xi (to plot)
ROUT1 = zeros(lpw, 1, nrows); %min loss of AC design
RTsigma = zeros(N, 1, nrows);
RESIGMA = zeros(N,lpw, nrows);
Rnuloss = zeros(Nl, lpw, nrows);
RnuSample = zeros(1, lpw, nrows);
RCxpoints = zeros(ds, lpw, nrows);

for row = 1:nrows
    [OUT0, OUT1, Tsigma, ESIGMA, nuSample, nuloss, Cxpoints] = adaptive(ds, plotwhich, theta0, eta, tau, SIGMA(row,:), Nl, N, X, initheta, ix);
    ROUT0(:,:,row) = OUT0; %Design
    ROUT1(:,:,row) = OUT1;
    RTsigma(:,:,row) = Tsigma;
    RESIGMA(:,:,row) = ESIGMA;
    Rnuloss(:,:,row) = nuloss;
    RnuSample(:,:,row) = nuSample + 2;
    RCxpoints(:,:,row) = Cxpoints;
end

%%Computing Real losses for Adaptive designs
for row = 1:nrows
    for count = 1:lpw
        Atrueloss(count, row) = round(maxamse(ROUT0(:,count,row), Realsigma(row,:)', TQ1, N, plotwhich(count)), 4);
    end
end

Maxt = 0;
for count = 1:lpw
       Maxc = max(RTsigma(:,1,count));
       Maxt = max(Maxc, Maxt);
end

Maxe = 0;
for row = 1:nrows
   for count = 1:lpw
       Maxcr = max(RESIGMA(:, count, row));
       Maxe = max(Maxcr, Maxe);
   end
end

%% ====================== Adaptive design points ========================
function [OUT0, OUT1, Tsigma, ESIGMA, nuSample, nuloss, Cxpoints] = adaptive(ds, plotwhich, theta0, eta, tau, sigma, Nl, N, X, initheta, ix)

lpw = length(plotwhich);
nuloss = zeros(Nl, lpw);
nuSample = zeros(1, lpw);
Cxpoints = zeros(ds, lpw);
OUT0 = zeros(N, lpw); % X versus xi (to plot)
OUT1 = zeros(lpw, 1); %min loss of AC design

%Data generation
m = 20; 
F0all = (theta0(1)*X')./(theta0(2)+X'); % Generated using theta0
%%Derivative of F w.r.t theta
fall = [X'./(theta0(2)+X'), -(theta0(1)*X')./((theta0(2)+X').^2)]; 

Xv = zeros(N,m);
Yv = zeros(N,m); %Generate repeated responses for population

for j = 1:m
 Xv(:,j) = X';
end 

%To store sigma at the final stage of estimation process
ESIGMA = zeros(N,lpw);

% Standardized Variance function
Tsigma = sqrt(sigma'/mean(sigma));

%Construction of model error delta
[Q, ~] = qr(fall);
p = size(fall,2);
Q2 = Q(:,(p+1):N);
E = normrnd(0,1,[N-p, 1]); % A random vector with size N-p from a standard normal distribution
c = E/norm(E); %%norm c equals to 1
deltax = eta*sqrt(N)*Q2*c;

% Generating errors from standard normal distribution
a = 2; %sigma epsilon
b = norminv(tau);

for j = 1:m
 errors = (a*Tsigma).*randn(N,1) - Tsigma*a*b;
 Yv(:,j) = F0all + errors;
end

Yp = Yv(:);

%AC: Adaptive design to construct design points until convergence
for count = 1:lpw
[xiAC, lossxiAC, Loss, nG, esigma, x] = adaptconvergence(ds, plotwhich(count), tau, N, X, Yv, initheta, F0all, deltax, Tsigma, a, b, Nl, Yp, ix);
OUT0(:, count) = xiAC; % A design
OUT1(count, 1) = lossxiAC; % Losses of Adaptive designs
nuloss(:, count) = Loss;
nuSample(1, count) = nG; 
ESIGMA(:,count) = esigma;
Cxpoints(:,count) = sort(x);
end

%% ============================= Adaptive design until convergence [AC] =============================
function [xi, lossxi, Loss, nG, esigma, x] = adaptconvergence(ds, nu, tau, N, X, Yv, initheta, F0all, deltax, Tsigma, a, b, Nl, Yp, ix)

eIND = (1:N)';
IX = zeros(N, 2);
IX(:,1) = eIND;
IX(:,2) = X';

h = 40; %Bandwidth size

Ind = 1:N;
Loss = zeros(Nl, 1);

x = ix;
n1 = length(x);
n = n1;
sInd = Ind(ismember(X',x)); % Initial sample indices
y = zeros(n,1);
for j = 1:n
    y(j) = datasample(Yv(sInd(j),:), 1);
end

fA = ones(n,1); 
freqA = zeros(N,1);
freqA(ismember(X,x))= fA;
xi = freqA/n; %initial design

ngeneration = 0;

while ngeneration < ds - n1
ngeneration = ngeneration + 1;

% fminsearch
% Quntile objective function that is required to minimize
fminfun = @(theta) sum(((y - theta(1)*x./(theta(2) + x))).*(tau - vindicator(y - theta(1)*x./(theta(2) + x))));
thetahat = fminsearch(fminfun, initheta); % estimates of theta

%Estimation of scale functions
esigma = estsigma(x, y, N, IX, h, thetahat);

% Estiamtes for Derivative of F w.r.t theta
efall = [X'./(thetahat(2)+X'), -(thetahat(1)*X')./((thetahat(2)+X').^2)]; 
p = size(efall,2);

[Q, ~] = qr(efall);
Q1 = Q(:,1:p);

[xinew, minlossnew, xnew, ynew] = adaptsequential(x, y, xi, esigma, Q1, N, nu, n, X, Yv, F0all, deltax, Tsigma, a, b, Yp);
    
    xi = xinew;
    x = xnew;
    y = ynew;
    n = n + 1;
    Loss(ngeneration) = minlossnew;
end
% required output
lossxi = minlossnew; 
nG = ngeneration;

%% =============      adaptsequential       =======================================================
%% ==============One design point selection at each time using previous information[adaptive]=======================
function [xinew, minlossnew, xnew, ynew] = adaptsequential(x, y, xi, esigma, Q1, N, nu, n, X, Yv, F0all, deltax, Tsigma, a, b, Yp)
    Index = 1:N;
    nextlosses = zeros(N,1);
    for i = 1:N
    fnext = zeros(N,1); 
    fnext(i) = 1; 
    xinext = (n*xi + fnext)./(n+1);
    nextlosses(i) = maxamse(xinext, esigma, Q1, N, nu); 
    end
    minlossnew = min(nextlosses); 
    minlossind = Index(ismember(nextlosses, minlossnew));
    istar = datasample(minlossind, 1);
    fnext = zeros(N,1); 
    fnext(istar) = 1; 
    xinew = (n*xi + fnext)./(n+1);
    xnew = [x', X(istar)]';
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 oneerror = ((a*Tsigma(istar)).*randn(1,1)) - (Tsigma(istar)*a*b);
 yadd = F0all(istar) + ((deltax(istar)/sqrt(n+1))*Tsigma(istar)) + oneerror;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if yadd > 0
           ynew = [y', yadd]'; 
    else
        pidx = Yp > 0;
        allpos = Yp(pidx);
        allmin = min(allpos);
        Yvistar = Yv(istar,:);
        idx = Yvistar > 0; % create logical index
        posistar = [Yvistar(idx), allmin];
        yadd = datasample(posistar, 1);
        ynew = [y', yadd]';
    end
   
%% Estimation of scales %%%%%%%%%%%%%%%%%%%%%%
function [esigma] = estsigma(x, y, N, IX, h, thetahat)

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

rx = [ones(nonz, 1) regdata(:,1)];
rcof = regress(regdata(:,2), rx);
ax = rcof(1);
bx = rcof(2);
mv = min(regdata(:,2))/10;

for j = 1:N
        sdx(j) = ax + (bx*IX(j, 2));
        if sdx(j) <= 0
            sdx(j) = mv;
        end
end

esigma = sdx/sqrt(mean(sdx.^2));

%% ============================= maxamse =============================
function L = maxamse(xi, sigma, Q1, N, nu)
p = size(Q1,2);
D = diag(xi);
Ip = eye(p);
xds = xi./sigma;
M = diag(xds);
Msq = diag(xds.^2);
V = (Q1'*D)*Q1;
H = (Q1'*M)*Q1;
r = rcond(H);
if r > 1e-15
    Hinv = H\Ip;
    A = (Hinv*V)*Hinv;
    L1 = trace(A);
    S = (Q1'*Msq)*Q1;
    B = real((Hinv*S)*Hinv);
    [~, L2] = eigs(B,1);
else
    L1 = 100;
    L2 = 100;
end
L = (1 - nu)*L1 + nu*L2;
L = L/N;

%% ============================= Scale =============================
function result = Scale(xv, x, h) 
    % xv is the input (1 x 1)
    % x is the sample design points (n x 1)
    % y is the sample responses (n x 1)
    % thetahat is the tau^th quantile estimates using sample (p x 1)
    % h is the bandwidth 
    kweight = exp((-0.5)*(((xv - x)./h).^2)); % Standard Gaussian kernel
    result = kweight./sum(kweight); % Standardized weights
    
%% ============================= vindicator =============================
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


%% ============================= Fun =============================
function out = Fun(x, theta)  
    % Forms the vector of nonlinear, Michaelis-Menten,  responses
    % theta is parameter (2 x 1)
    % x is a design vector (n x 1), where n is the sample size
    % out (n x 1) 
out = theta(1)*x./(theta(2) + x);

%% ============================= closest =============================
function xmin = closest(X,y)
% finds the closest member of X to y.
dist = abs(X-y);
d0 = min(dist);
xmin = X(dist==d0);
if length(xmin) > 1
    xmin = xmin(1);
end
