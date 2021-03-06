function theta = getStandardParameters(theta,type,widthalpha)
% this function transforms a parameter given in threshold, width to a
% standard parameterization for comparison purposes
%function theta = getStandardParameters(theta,type,[widthalpha])
% theta = parameter
% type  = Name of the Sigmoid
%
% if you changed the widthalpha you should pass it as an additional
% argument
%
% norm/gauss/normal
%       theta(1) = threshold -> threshold
%       theta(2) = width     -> standard deviation
% logistic
%       theta(1) = threshold -> threshold
%       theta(2) = width     -> slope at threshold
% Weibull/weibull
%       theta(1) = threshold -> scale 
%       theta(2) = width     -> shape parameter
% gumbel & rgumbel distributions
%       theta(1) = threshold -> mode
%       theta(2) = width     -> scale
% tdist
%       theta(1) = threshold -> threshold/mean
%       theta(2) = width     -> standard deviation

if ~exist('widthalpha','var'), widthalpha=.05; end

if isstruct(theta)
    widthalpha = theta.options.widthalpha;
    type  = theta.options.sigmoidName;
    if theta.options.threshPC~=.5
        if theta.options.logspace
            theta.Fit(1) = log(getThreshold(theta, .5, true));
        else
            theta.Fit(1) = getThreshold(theta, .5, true);
        end
    end
    theta = theta.Fit;
end

assert(logical(exist('type','var')),'You need to specify which kind of sigmoid you fit');

switch type
    case {'norm','gauss','normal'}
        theta(1) = theta(1);
        C        = my_norminv(1-widthalpha,0,1) - my_norminv(widthalpha,0,1);
        theta(2) = theta(2)./C;
    case 'logistic'
        theta(1) = theta(1);
        theta(2) = 2*log(1./widthalpha-1)./theta(2);
    case {'Weibull','weibull'}
        C        = log(-log(widthalpha)) - log(-log(1-widthalpha));
        shape    = C./theta(2);
        scale    = exp(C./ theta(2) .* (-theta(1)) + log(-log(.5)));
        scale    = exp(log(1/scale)/shape); %Wikipediascale
        theta(1) = scale;
        theta(2) = shape;
    case {'gumbel'}
        % note that gumbel and reversed gumbel definitions are sometimes swapped
        % and sometimes called extreme value distributions
        C        = log(-log(widthalpha)) - log(-log(1-widthalpha));
        theta(2) = theta(2)/C;
        theta(1) = theta(1)-theta(2).*log(-log(.5));
    case {'rgumbel'}
        C      = log(-log(1-widthalpha)) - log(-log(widthalpha));
        theta(2) = -theta(2)./C;
        theta(1) = theta(1)+theta(2).*log(-log(.5));
    case {'tdist','student','heavytail'} 
        C      = (my_t1icdf(1-widthalpha) - my_t1icdf(widthalpha));    
        theta(1) = theta(1);
        theta(2) = theta(2)./C;
        
end