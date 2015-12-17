function [powr, lambda] = cor_power(x,alpha,lambda,w)
% Based on implementation by Miika Ahdesmaki and Korbinian Strimmer in R
% package sda (2013-11-19).
%
% Jasper Engel 17-12-2015.

[n,p] = size(x);

if nargin < 4 || isempty(w)
    w = ones(n,1)*(1/n);% Can be adjusted later if I want to put in weighting for samples.
end

if nargin < 3 || isempty(lambda)
    lambda = estimate_lambda(x,w);
elseif lambda < 0 
    lambda = 0;
elseif lambda > 1
    lambda = 1;
end

% Autoscale
mx = mean(x);
sx = x - ones(n, 1)*mx;
sd = ones(n, 1)*std(sx);
xs = sx./sd; 

% Bias correction factor
h1 = 1/(1-sum(w.*w)); % For w = 1/n this equals the usual h1 = n/(n-1)

if lambda == 1 && alpha == 0% Result in both cases is the identity matrix
    powr = eye(p);
elseif alpha == 1 % No SVD in this case
    r0 = h1*xs'*diag(w)*xs;
    powr = (1-lambda)*r0;
    powr(logical(eye(p))) = 1;
else
    zeros = var(xs) == 0;
    
    [svdxsu,svdxss,svdxsv] = svdecon(xs,1e-7);
    
    m = length(diag(svdxss)); % Rank of xs
    
    UTWU = svdxsu'*diag(w)*svdxsu;
    C = svdxss*UTWU*svdxss;
    C = (1-lambda)*h1*C;
    
    C = (C+C')/2; % Symmetrize for numerical reasons
    
    if lambda == 0 % Use eigenvalue decomposition computing the matrix power
        powr = svdxsv*mpower(C,alpha)*svdxsv';
    else 
        F = eye(m)-mpower(C/lambda+eye(m),alpha);
        powr = (eye(p) - svdxsv*F*svdxsv').*(lambda)^alpha;
    end
    
    % Set all diagonal entries corresponding to zero-variance variables to
    % 1 
    powr(diag(zeros)) = 1;
end