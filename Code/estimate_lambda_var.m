function lambda_var = estimate_lambda_var(x,w)
% Based on implementation by Miika Ahdesmaki and Korbinian Strimmer in R
% package sda (2013-11-19).
%
% Jasper Engel 17-12-2015.

[n,p] = size(x);

if nargin < 2
    w = ones(n,1)*(1/n);% Can be adjusted later if I want to put in weighting for samples.
end

% Bias correction factors 
w2 = sum(w.*w); % for w = 1/n this equals 1/n
h1 = 1/(1-w2); % for w = 1/n this equals the usual h1 = n/(n-1)
h1w2 = w2/(1-w2); % for w = 1/n this equals 1/(n-1)

% Center the data
m = mean(x);
xc = x - ones(n, 1)*m;

% Compute empirical variances
v = h1.*sum((w*ones(1, p)).*xc.^2);

% Compute shrinkage target
target = median(v);

zz = xc.^2;
q1 = sum(zz.*(w*ones(1, p)));
q2 = sum((zz.^2).*(w*ones(1, p)))-q1.^2;
numerator = sum(q2);
denominator = sum((q1 - target/h1).^2);

if denominator == 0
    lambda_var = 1;
else
    lambda_var = min(1,numerator/denominator*h1w2);
end