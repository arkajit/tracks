% ----------------------------------------------
% backward propagation
% Generalized to support HMMS with continuous emissions

function [log_b] = backward(x,hmm) 
    
k = size(hmm.log_t, 1); % number of states
m = length(x); % length of sequence

log_b = zeros(k, m); 
log_b(:,m) = 0; % b(:,m) = 1

for j=m-1:-1:1 	% for each timestep
	log_e = zeros(k, 1);
	for s=1:k
		log_e(s) = hmm.log_E(s, x(j+1));
	end
	log_E = repmat(log_e, 1, k);
	% log_E is a (k x k) matrix; each column is the same and
	% equal to hmm.log_E(:,x(j+1));

	log_B = repmat(log_b(:,j+1), 1, k); 
	% log_B is a (k x k) matrix; each column is the same and
	% equal to log_b(:,j+1);

	log_b(:,j) = log_sum_exp(hmm.log_T + log_E' + log_B', 2);
end;
