% ---------------------------------
% turns log-probabilities into probabilities 

function [p] = norm_exp(log_p) 

p = exp(log_p(:)-max(log_p(:)));  % k*m-length col vector
p = reshape(p./sum(p), size(log_p) ); % normalize and reshape into kxm matrix

