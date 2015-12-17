function [Cov, lambda, lambda_var] = cov_power(x,alpha,lambda_cor,lambda_var,w)
% Based on implementation by Miika Ahdesmaki and Korbinian Strimmer in R
% package sda (2013-11-19).
%
% Jasper Engel 17-12-2015.

[n,~] = size(x);

if nargin < 5 || isempty(w)
    w = ones(n,1)*(1/n);% Can be adjusted later if I want to put in weighting for samples.
end

if nargin < 4 || isempty(lambda_var)
    lambda_var = estimate_lambda_var(x,w);
elseif lambda_var < 0 
    lambda_var = 0;
elseif lambda_var > 1
    lambda_var = 1;
end

if nargin < 3 || isempty(lambda_cor)
    lambda = estimate_lambda(x,w);
else
    lambda = lambda_cor;
    if lambda < 0
        lambda = 0;
    elseif lambda > 1
        lambda = 1;
    end
end

% Shrinkage scale factors
[vs, ~] = svar(x,lambda_var, w);
sc = sqrt(vs);

% Shrinkage correlation
[powr_cor, ~] = cor_power(x,alpha,lambda,w);

% Shrinkage covariance 
Cov = (sc'*sc).*powr_cor;