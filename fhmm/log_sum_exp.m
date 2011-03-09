% ----------------------------------------------
% sums probabilities stored as log-probabilities,
% analogously to MATLAB sum() function.
% Returns the answer as a log-probability

function [log_p] = log_sum_exp(log_P,ind)

log_max  = max(log_P,[],ind); 
log_Pmax = repmat(log_max,size(log_P)./size(log_max));
I = find(log_Pmax>-Inf); % shouldn't subtract -Inf 
log_Pmax(I) = log_P(I) - log_Pmax(I); 
log_p = log_max + log(sum(exp(log_Pmax),ind));

