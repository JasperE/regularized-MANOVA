function [vs, lambda_var] = svar(x,lambda_var, w)
% Shrinkage of variances. 
% Based on implementation by Miika Ahdesmaki and Korbinian Strimmer in R
% package sda (2013-11-19).
%
% Ref: 
%
% Jasper Engel 17-12-2015.

[n,p] = size(x);

if nargin < 3 || isempty(w)
    w = ones(n,1)*(1/n);% Can be adjusted later if I want to put in weighting for samples.
end

if nargin < 2 || isempty(lambda_var)
    lambda_var = estimate_lambda_var(x,w);
elseif lambda_var < 0 
    lambda_var = 0;
elseif lambda_var > 1
    lambda_var = 1;
end

% Bias correction factor % for w = 1/n this equals the usual h1 = n/(n-1)
h1 = 1/(1-sum(w.*w));

% Center the data
m = mean(x);
xc = x - ones(n, 1)*m;

% Compute empirical variances
v = h1.*sum((w*ones(1, p)).*xc.^2);

% Compute shrinkage target
target = median(v);

% Shrinkage estimate 
vs = lambda_var*target + (1-lambda_var)*v;