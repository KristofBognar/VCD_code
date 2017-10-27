function [w_mean, w_sigma, w_sigma_corrected] = get_weighted_mean(x, sigma)
% Written Cristen Adams
% February 2011
% 
% This function calculates the weighted mean of the given values and errors
%
%

L = length(x);
if length(sigma) ~= length(x),
    error('x and sigma must be vectors of the same length')
end

w_mean = sum(x ./ (sigma.^2)) / sum( 1 ./ (sigma.^2));
w_sigma = sqrt(1 ./ sum( 1 ./ (sigma.^2))   );
w_sigma_corrected = sqrt(1 ./ sum( 1 ./ (sigma.^2)) *...
    1/(L-1) * sum( (x - mean(x)).^2 ./ (sigma.^2) ) );