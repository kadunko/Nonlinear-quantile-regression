%% ====================== Main function to choose design points by SL method ========================
function [SLTrueSigma1Table, SLTrueSigma2Table, SLTrueSigma3Table,... 
    SLTrueSigma4Table, FTrueSigma1Table, FTrueSigma2Table,...
    FTrueSigma3Table, FTrueSigma4Table] = nonlinearsequential(ds, plotwhich, n0)
% ds is the desired sample size.
% plotwhich is the desired values of nu
% n0 is the desired number of fixed points
% e.g: [SLTrueSigma1Table, SLTrueSigma2Table, SLTrueSigma3Table, SLTrueSigma4Table, FTrueSigma1Table, FTrueSigma2Table, FTrueSigma3Table, FTrueSigma4Table] = nonlinearsequential(1000, [.2 .35 .65 .95], 30)

format short g
close ALL
rng(200); % For reproducibility

atheta0 = [57.98, 46.43]'; % Assumed initial estimates for SL design

XL = 5; % Lower limit of design space
XU = 400; % Upper limit of design space
N = 100; % Total number of points in design space
X = linspace(XL, XU, N); % Design space

%%Derivative of F w.r.t theta for SL design
fall = [X'./(atheta0(2)+X'), -(atheta0(1)*X')./((atheta0(2)+X').^2)]; 

Nl = 5000; %Preallocate size to save losses
p = size(fall,2); %The number of model parameters

[Q, ~] = qr(fall); % QR decomposition
Q1 = Q(:,1:p);

% Real sigma functions
Rvar1 = ones(1,N);
Rvar2 = (X + X(N)).^2;
Rvar3 = X + X(N);
Rvar4 = ((2*X) + (0.5*X(N))).^2;

% Standardized sigma functions for satisfying the condition (2ii)
Rsigma1 = sqrt(Rvar1/mean(Rvar1)); 
Rsigma2 = sqrt(Rvar2/mean(Rvar2));
Rsigma3 = sqrt(Rvar3/mean(Rvar3));
Rsigma4 = sqrt(Rvar4/mean(Rvar4));

Rsigma = [Rsigma1; Rsigma2; Rsigma3; Rsigma4];

% Assumed Variance functions to construct SL designs
var1 = ones(1,N);
var2 = (X + X(N)).^2;
var3 = X + X(N);
var4 = ((2*X) + (0.5*X(N))).^2;

VAR = [var1; var2; var3; var4];
nrows = size(VAR, 1);
lpw = length(plotwhich);
RCFixpoints = zeros(n0, lpw, nrows); % n0-point Fixed design

SLtrueloss = zeros(nrows, lpw, nrows);
Ftrueloss = zeros(nrows, lpw, nrows);

LossSL = zeros(lpw, 2, nrows); %minlosses of SL designs with nu
LossFix = zeros(lpw, 2, nrows); %minlosses of n0-point designs with nu
Mloss = zeros(Nl, lpw, nrows);
Msample = zeros(nrows, lpw);
ROUT0 = zeros(N, lpw, nrows); % To save designs
RCFdesign = zeros(N, lpw, nrows); % n0 points for fixed designs

for row = 1:nrows
    [OUT0, OUT1, OUT2, nuloss, nuSample, CFdesign, CFixpoints] = nonlinear(X, Q1, ds, N, n0,  VAR(row,:), plotwhich, nrows, row, Nl);
    ROUT0(:,:,row) = OUT0;
    LossSL(:,:,row) = OUT1;
    LossFix(:,:,row) = OUT2;
    Mloss(:,:,row) = nuloss;
    Msample(row,:) = nuSample;
    RCFdesign(:, :, row) = CFdesign;
    RCFixpoints(:, :, row) = CFixpoints;
end

%%Real losses for SL designs
for nrow = 1:nrows
    for row = 1:nrows
        for count = 1:lpw
        SLtrueloss(row, count, nrow) = round(maxamse(ROUT0(:,count,row), Rsigma(nrow,:), Q1, N, plotwhich(count)), 4);
        end
    end
end

% SL Loss Tables
rowNames = {'nu = 0.20', 'nu = 0.35', 'nu = 0.65', 'nu = 0.95'};
colNames = {'SL1', 'SL2', 'SL3', 'SL4'};
SLTrueSigma1Table = array2table(SLtrueloss(:, :, 1)','RowNames',rowNames,'VariableNames',colNames);
SLTrueSigma2Table = array2table(SLtrueloss(:, :, 2)','RowNames',rowNames,'VariableNames',colNames);
SLTrueSigma3Table = array2table(SLtrueloss(:, :, 3)','RowNames',rowNames,'VariableNames',colNames);
SLTrueSigma4Table = array2table(SLtrueloss(:, :, 4)','RowNames',rowNames,'VariableNames',colNames);

%%Real losses for Fixed n0-point designs
for nrow = 1:nrows
    for row = 1:nrows
        for count = 1:lpw
        Ftrueloss(row, count, nrow) = round(maxamse(RCFdesign(:,count,row), Rsigma(nrow,:), Q1, N, plotwhich(count)), 4);
        end
    end
end

FTrueSigma1Table = array2table(Ftrueloss(:, :, 1)','RowNames',rowNames,'VariableNames',colNames);
FTrueSigma2Table = array2table(Ftrueloss(:, :, 2)','RowNames',rowNames,'VariableNames',colNames);
FTrueSigma3Table = array2table(Ftrueloss(:, :, 3)','RowNames',rowNames,'VariableNames',colNames);
FTrueSigma4Table = array2table(Ftrueloss(:, :, 4)','RowNames',rowNames,'VariableNames',colNames);

%% ========= nonlinear =====================================
function [OUT0, OUT1, OUT2, nuloss, nuSample, CFdesign, CFixpoints] = nonlinear(X, Q1, ds, N, n0, var, plotwhich, nrows, row, Nl)

lpw = length(plotwhich);
OUT0 = zeros(N, lpw); % To save designs
OUT1 = zeros(lpw, 2); %minlosses of SL designs
OUT2 = zeros(lpw, 2); %minlosses of n-point designs
CFdesign = zeros(N, lpw); % n0-point Fixed design
CFixpoints = zeros(n0, lpw); % n0-point Fixed design

nuloss = zeros(Nl, lpw);
nuSample = zeros(1, lpw);

sigma = sqrt(var/mean(var)); % mean(sigma.^2) = 1

for count = 1:lpw
[xiSL, LossSL, Loss, nSL] = sequential(X, sigma, Q1, N, ds, plotwhich(count), Nl);
OUT0( :, count) = xiSL;
OUT1(count, 2) = LossSL;
OUT1(count, 1) = plotwhich(count);
nuloss(:,count) = Loss;
nuSample(count) = nSL;
end

figure(1) %% Plots of designs
for count = 1:length(plotwhich)
    xi = zeros(N,1);
    Fixpoints = zeros(1,n0);
    subplot(nrows, lpw, lpw*(row - 1) + count);
    xistar = OUT0( :, count);
    mm1 = 1.2*max(xistar);
    mm2 = 1.5*max(xistar);
    mm = min(mm1,mm2);

    cdf = cumsum(xistar);
    bar(X, xistar, 'EdgeColor', 'b', 'FaceColor', 'none');
    axis([1 400 -.05*mm mm1]);
    hold on
    n = 0;
    for i = 1:n0
        vec = abs((i-.5)/n0 - cdf);
        [~,I] = min(vec);
        Fixpoints(i) = X(I);
        fnext = zeros(N,1); 
        fnext(I) = 1; 
        xi = (n*xi + fnext)/(n+1);
        n = n + 1;
    end
    CFdesign(:,count) = xi;
    
%Computing loss of n-point design
    L = maxamse(xi, sigma, Q1, N, plotwhich(count));
    OUT2(count, 2) = L;
    OUT2(count, 1) = plotwhich(count);
    
	scatter(Fixpoints, zeros(1,n0) - .05*mm, 10, 'filled', 'MarkerFaceColor', 'r');
    hold off
    
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)'];
    
    group = 4*(row - 1) + count;
    xlabel(labels((3*(group-1)+1):(3*(group-1)+3)));
	ylabel('\xi(x)');
   CFixpoints(:,count) = Fixpoints';
end

figure(2) %% Plots of designs
for count = 1:length(plotwhich)
    subplot(nrows, lpw, lpw*(row - 1) + count);
    plot(3:nuSample(count), nuloss(1:nuSample(count)-2,count));
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)'];
    
    group = 4*(row - 1) + count;
    xlabel(labels((3*(group-1)+1):(3*(group-1)+3)));
	ylabel('Loss');
end

%% ============================= sequential local design until convergence [SL]==========================
function [xiSL, LossSL, Loss, nSL] = sequential(X, sigma, Q1, N, ds, nu, Nl)
%%Initial design
n = 2; % initial design size
XD = [X(1), X(N)]';
fD = ones(n,1); 
freqD = zeros(N,1);
freqD(ismember(X,XD))= fD;
xi = freqD/n;

Loss = zeros(Nl,1);
Index = 1:N;

while n < ds
    nextlosses = zeros(N,1);
    for i = 1:N
    fnext = zeros(N,1); 
    fnext(i) = 1; 
    xinext = (n*xi + fnext)/(n+1);
    nextlosses(i) = maxamse(xinext, sigma, Q1, N, nu); 
    end
    minlossnew = min(nextlosses); 
    minlossind = Index(ismember(nextlosses, minlossnew));
    istar = datasample(minlossind, 1);
   
    fnext = zeros(N,1); 
    fnext(istar) = 1; 
    xi = (n*xi + fnext)/(n+1);
    n = n + 1;
    Loss(n-2) = minlossnew;
end

%Final output for sequential design
xiSL = xi;
L = maxamse(xi, sigma, Q1, N, nu);
LossSL = L; 
nSL = n;

%% ============================= maxamse =============================
function L = maxamse(xi, sigma, Q1, N, nu)
p = size(Q1,2);
D = diag(xi);
Ip = eye(p);
xds = xi./sigma';
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