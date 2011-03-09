
% ----------------------------------------------
% forward propagation 

function [log_a] = forward(x,hmm) 
    
k = size(hmm.log_t, 1); % number of states
m = length(x);         % length of sequence

log_a = zeros(k, m); 
for s=1:k 	% for each state
	log_a(s, 1) = hmm.log_t(s) + hmm.log_E(s, x(1)); 
end;

for j=2:m		% for each timestep
	log_A = repmat(log_a(:,j-1), 1, k); 
	%log_A is a (k x k) matrix where each column is the same and 
	% equal to log_a(:,j-1); 
	prev = log_sum_exp(log_A + hmm.log_T, 1)';
	for s=1:k 	% for each state
		log_a(s,j) = prev(s) + hmm.log_E(s, x(j));
	end;
end;
