% ---------------------------------
% turns log-probabilities into probabilities 

function [p] = norm_exp(log_p) 

p = exp(log_p(:)-max(log_p(:))); 
p = reshape(p./sum(p), size(log_p) );

