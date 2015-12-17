function stats = rmanova(y,options,varargin)
%% reguarlized Multivariate Analysis of Variance (rMANOVA)
%
% Input 
%   y: data matrix (samples x variables)
%   options: options structure
%   varargin: labels for each factor of interest (one column-vector per
%   factor)
%
% Output
%   stats: rmanova(y,options,varargin) returns the following fields:
%       stats.data: results fit of data by Eq. 2 from rMANOVA paper
%       stats.info: information regarding factors used in model; d.f.
%       correspond to d.f. of classic MANOVA Model
%       stats.dmat: design matrix that was used
%       stats.shrink: shrinkage parameter used by model
%       stats.test: Estimated test statistics for the factors and
%       interactions of interest. Currently Wilk's Lambda, Pillai's Trace,
%       Hotelling-Lawley Trace and Roy's Max Root are supported. The first
%       index corresponds to the tests and the second index to factors and
%       interactions (specified in stats.info).
%       stats.canon_val.std: Standardized canonical variate found for the
%       different tests. First and second index refer to the variables and
%       CV's, respectively. Third index is main effect or interaction
%       (specified in stats.info).
%       stats.scores: Data projected onto CV's. Can be used to construct a
%       scores plot. First index is sample, second index is CV, third index
%       is main effect or interaction. Note that for high dimensional data
%       the score plot typically presents a too optimistic view of the
%       separation between the conditions for a specific factor. In that
%       case it might be useful to project the data onto e.g. the top 250
%       variables specified by the CV.     
%
% The default options structure can be loaded via the command
%   options: rmanova('options')
%
% The options structure includes the following fields:
%   vartype: indicates what type of coding of the design matrix should be
%   used. Default is 1: sum-to-zero design. (nominal encoding) vartype 0
%   van be used to indicate a continuous variable (i.e. MANCOVA)
%   If var_type is present it must be the same lenght as varargin or be a
%   scalar. If it is a scalar it applies to all variables.
%
%   model: The model is either a scalar or a matrix of terms. If it is a
%   scalar it refers to the degree of the model. 1 indicates a 1st degree
%   model which contains only main effects. A 2nd degree model contains
%   main effects and all pairwise interactions. For more complex models a
%   matrix of terms can be used. Each row of the matrix represents one term
%   in the model. Each column represents a predictor variable. A column
%   containing a non-zero value will be included in the given term of the
%   model. A row with all zeros is a global intercept term. If present it
%   must be the first row.
%
%   sstype: type of sums-of-squares that should be computed. Default is
%   type III.
% 
%   lambda: shrinkage parameter for correlation matrix. lambda = []
%   indicates that a James-Stein estimate of the within-group scatter
%   matrix should be estimated. Fixed values between 0 and 1 can also be
%   used as input. lambda = 1 corresponds to the ASCA model.
%
%   target: shrinkage target. The following targets can be used:
%       'average_var': diagonal with average variance
%       'unique_var': diagonal with unique variance
%
% Dependencies:
%
% glm: function from linstats toolbox 2006b by M. Boedigheimer
% (http://uk.mathworks.com/matlabcentral/fileexchange/13493-linstats-2006b)
% to specify a generalized linear model
%
% See examples_rmanova.m for an example
%
% References: Engel, J., et al. "Regularized MANOVA (rMANOVA) in untargeted
% metabolomics." Analytica chimica acta 899 (2015): 1-12.
%
% Jasper Engel 17-12-2015.

%% Code
if nargin == 1 % return options structure
    load options
    stats = options;  
else
    if isempty(options)
        load options
    end
    
    %% A: define linear model and design matrix
    %1. Create designmatrix using function from the linstats_2006b toolbox
    glm = encode(y, options.vartype, options.model, varargin{:});
    for i = 1:length(varargin) %2. save factornames
        glm.var_names{i,2} = inputname(2+i);
    end
    stats.info.factors = glm.var_names;
    stats.info.levels = glm.level_names;
    stats.info.model = glm.coeff_names;
    %3. Save labels main effects and interactions  
    for i = 1:length(unique(glm.terms))-1;
        stats.info.labels{i}=glm.dmat(:,glm.terms==i);
    end
    
    %% B: fit linear model
    stats_temp = fit_lm(glm, options);
    clear glm; 
    stats.data.raw = stats_temp.y;
    stats.data.fit = stats_temp.yhat; % unique(stats.data.fit(:,:,i)) gives the estimated means for the ith term in the model. Note that the means may seem slightly different due to rounding errors. 
    stats.data.res = stats_temp.resid;
    stats.dmat = stats_temp.dmat;
    stats.nmodels = stats_temp.nmodels;
    stats.ntests = stats_temp.ntests;
    
    % Degrees of freedom in normal MANOVA model; only accuracte for balanced data
    stats.info.df.total = stats_temp.dft; %Total degrees of freedom
    stats.info.df.reg = stats_temp.dfr; %Regresion degrees of freedom
    stats.info.df.res = stats_temp.dfe; %Residual degrees of freedom
    stats.info.df.hyp = stats_temp.dfh; %Hypothesis degrees of freedom 
    
    %% C: Estimate shrinkage parameters
    switch options.target
        case{'average_var'}
            stats.shrink.lambda_var = 1;
        case{'unique_var'}
            stats.shrink.lambda_var = 0;
        case{'shrink_var'}
            stats.shrink.lambda_var = estimate_lambda_var(stats.data.res); 
    end
    
    if isempty(options.lambda) % Estimate correlation shrinkage
        stats.shrink.lambda = estimate_lambda(stats.data.res);
    else
        stats.shrink.lambda = options.lambda;
        if stats.shrink.lambda > 1
            stats.shrink.lambda = 1;
        elseif stats.shrink.lambda < 0
            stats.shrink.lambda = 0;
        end
    end
    
    if options.extra == 1 % Dont save canonical variates and their scores (useful during permutation testing).
        [Cov_res, ~, ~] = cov_power(stats.data.res,1,stats.shrink.lambda,stats.shrink.lambda_var);
        Cov_res = Cov_res.*stats.info.df.total./stats.info.df.res;
        vs = diag(sqrt(Cov_res));
        Cor_res = Cov_res./((vs*vs'));
        clear vs;
    end
    
    %% D: Calculate eigenvectors and eigenvalues
    for i = 1:stats.ntests
        [eigenval, eigenvec, stdval] = invwb_prod(stats.data.res,squeeze(stats.data.fit(:,:,i)),stats.shrink.lambda,stats.shrink.lambda_var);
        stats.eigval(1:length(eigenval),i) = eigenval;
        
        if options.extra == 1
            stdval = stdval./sqrt(stats.info.df.res); %Equal to sqrt(diag(Cov_res))
            
            % Scale standardized eigenvectors so that v*Sigma*v = 1, where
            % sigma is W*/(n-1);
            % Not exactly the covariance matrix since total degrees of freedom used.
            vs = diag((eigenvec' * Cor_res * eigenvec))'; % Cov_res is correlation matrix in this case.
            vs(vs<=0) = 1;
            eigenvec = eigenvec ./ (ones(size(eigenvec,1), 1)*sqrt(vs));
            % Flip sign so that the average element of eigenvectors is positive
            k = (sum(eigenvec) < 0);
            eigenvec(:,k) = -eigenvec(:,k);
            stats.canon_var.std(:,1:size(eigenvec,2),i) = eigenvec;
            
            % Scale unstandardized eigenvectors to match output matlab manova1
            % function when no shrinkage is being applied.
            eigenvec_raw = eigenvec./(stdval(1,:)'*ones(1,size(eigenvec,2)));
            vs = diag((eigenvec_raw' * Cov_res * eigenvec_raw))';
            eigenvec_raw = eigenvec_raw ./(ones(size(eigenvec_raw,1), 1)*sqrt(vs));
            k = (sum(eigenvec_raw) < 0);
            eigenvec_raw(:,k) = -eigenvec_raw(:,k);
            stats.canon_var.raw(:,1:size(eigenvec_raw,2),i) = eigenvec_raw;
            
            % Project data onto eigenvectors and store scores
            % Scores of data projected onto canonical variates
            %y_center = (stats.data.raw - ones(size(stats.data.raw,1), 1)*mean(stats.data.raw))./stdval; %Should correct for other factors
            data_raw = squeeze(stats.data.fit(:,:,i)) + stats.data.res;
            y_center = (data_raw - ones(size(data_raw,1), 1)*mean(data_raw));
            y_std = y_center./stdval;
            stats.scores.raw(:,1:size(eigenvec_raw,2),i) = y_center*eigenvec_raw;
            stats.scores.std(:,1:size(eigenvec,2),i) = y_std*eigenvec;
        end
    end
    
    %% E: Calculate test statistics
    % Compuate test statistics all methods are functions of the eigenvalues
    % ('Wilk''s Lambda', 'Pillai''s Trace', 'Hotelling-Lawley Trace', 'Roy''s Max Root')
    for i = 1:stats.ntests
        e = stats.eigval(:,i); % Select eigenvalues of relevant test
        stats.tests(:,i)=[ prod(1./(1+e));  sum(e./(1+e)); sum(e);        max(e) ];
    end
    
    %% F: store additional output
    if options.extra == 1 % Save additional output       
        % Matrices of sums of squares and cross products
        stats.ss.total = (cov(stats.data.raw).*stats.info.df.total);
        stats.ss.res = Cov_res.*stats.info.df.res;
        stats.ss.hyp = zeros(size(stats.data.raw,2),size(stats.data.raw,2),stats.ntests);
        for i = 1:stats.ntests
            stats.ss.hyp(:,:,i) = stats.data.fit(:,:,i)' *stats.data.fit(:,:,i);
        end
    end
    clear stats_temp; clear powr_cov; clear powr_cor; clear Cov_res; clear Cor_res;
end

function [fit_stats] = fit_lm( glm, options)
% Fit linear model specified by glm using the sums of squares type
% indicated by sstype. 
sstype = options.sstype;
% -squares solution
ls_full =  lsr(glm.dmat,glm.y);
fit_stats.y = glm.y;
fit_stats.dmat = glm.dmat;

[n,q] = size(glm.y);

%% Between Group Sum Squares and Cross Products
% ss type I, II or III tests
[tests, reference] = gettests( glm, sstype);

nmodels = size(tests,1);
ntests  = length(reference);

fits = zeros(n,q,nmodels+1);

% leave out each term and calculate its marginal
% reduction in the sum of squared errors
for i = 1:nmodels
    j = ~ismember( glm.terms, find(tests(i,:) ));
    
    % testing the intercept is equivalent to comparing
    % two models: 1) with all terms except intercept to the
    % full model which includes all terms and the intercept
    ls =  lsr( glm.dmat(:,j), glm.y);
    fits(:,:,i) = ls.yhat; % reduced model
    dfh(i) = length(find(j == 0));
end;

% ssr for full model
% this isn't flexible in that it assumes nmodels+1 is full model
fits(:,:,nmodels+1) = ls_full.yhat;
yhat = fits(:,:,reference) - fits(:,:,1:ntests);

% save output
fit_stats.yhat = yhat(:,:,1:ntests); %use these to calculate the SSQ matrices
fit_stats.resid = ls_full.resid;     % use to calculate residual matrix
fit_stats.dft = ls_full.dft;
fit_stats.dfr = ls_full.dfr;
fit_stats.dfe = ls_full.dfe;
fit_stats.dfh = dfh(1:ntests); % ONLY correct when the data is balanced
fit_stats.nmodels = nmodels;
fit_stats.ntests = ntests;

function [ls] = lsr( dmat, y)
% Least squares regression fit.
ls.yhat = dmat*pinv(dmat'*dmat)*dmat'*y;
ls.resid = y - ls.yhat;
ls.dft = size(y,1) - 1;
ls.dfr = size(dmat,2)-1;
ls.dfe = ls.dft - ls.dfr;

