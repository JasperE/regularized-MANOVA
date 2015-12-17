function lambda = estimate_lambda(x,w)
% Based on implementation by Miika Ahdesmaki and Korbinian Strimmer in R
% package sda (2013-11-19).
%
% Jasper Engel 17-12-2015.


[n, p] = size(x);

if nargin < 2
    w = ones(n,1)*(1/n);% Can be adjusted later if I want to put in weighting for samples.
end

if p == 1
    lambda = 1;
end

% Autoscale
m = mean(x);
sx = x - ones(n, 1)*m;
sd = ones(n, 1)*std(sx);
xs = sx./sd; 

% Bias correction factors 
w2 = sum(w.*w); % for w = 1/n this equals 1/n
h1w2 = w2/(1-w2); % for w = 1/n this equals 1/(n-1)

sw = sqrt(w);

% Compute off-diagonal sums much more efficiently for n << p
xsw = xs.*(sw*ones(1, p));
[xswu,xsws,xswv] = svdecon(xsw,1e-7); % Adjust svdecon to only retain eigenvectors with positive eigenvalues
sE2R = sum(sum(xsw.*((xswu*xsws.^3)*xswv')))-sum(sum(xsw.^2).^2);
clear xsw; clear xswu; clear xsws; clear xswv; % Free memory
xs2w = (xs.^2).*(sw*ones(1, p));
sER2 = 2*sum(sum(xs2w(:,p-1:-1:1).*(cumsum(xs2w(:,p:-1:2),2))));
clear xs2w; % Free memory

denominator = sE2R;
numerator = sER2 - sE2R;

if denominator == 0
    lambda = 1;
else
    lambda = min(1,numerator/denominator *h1w2);
end