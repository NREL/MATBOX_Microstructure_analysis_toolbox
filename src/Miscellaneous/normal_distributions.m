function [x, pdf_, cdf_] = normal_distributions(equiprobability, mean_, sigma_, maxdeviation)
% [x, pdf_, cdf_] = normal_distributions(equiprobability, mean_, sigma_, maxdeviation)
% If equiprobability is true:
%   output is a flat distribution from mean_ - maxdeviation to mean_ + maxdeviation
% Otherwise:
%   output is a probability density function of a gaussian distribution with mean_ and sigma_
% Outputs:
% - x    : values
% - pdf_ : probability density function of x
% - cdf_ : cumulative distribution function of x

x = []; pdf_ = []; cdf_ =[];
% Univariate Gaussian Distribution
% Probability density function of a gaussian distribution
pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

if ~equiprobability && maxdeviation~=0
    x=linspace(mean_-maxdeviation, mean_+maxdeviation, 1000);
    pdf_=pdf_GD(x,mean_,sigma_);
    cdf_ = pdf2cdf(x,pdf_);
end

    function c = pdf2cdf(x,pdf)
        n=length(x);
        c=zeros(1,n);
        for i=2:1:n
            c(i) = trapz(x(1:i),pdf(1:i));
        end
        if c(end)<1
            c(end)=1;
        end
    end
end



