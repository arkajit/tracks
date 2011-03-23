function [hmm] = em_fhmm(X, hmm)
% em_fhmm - Expectation Maximization for a continuous Gaussian HMM
%		X - cell array of observed emissions, each vector one training example
%		hmm - initial guess for the parameters
%			requires t, T, means, stddevs fields to be set to init values
%			returns the MLE of hmm parameters including t, T, means, stddevs

[hmm, L] = one_em_iter(X, hmm);

L0 = L-abs(L);
while (L-L0 > abs(L)*1e-6)
	L0 = L;
	[hmm, L] = one_em_iter(X, hmm);
end

hmm.loglik = L;

% -----------------------------------------
function [hmm, loglik] = one_em_iter(X, hmm)
	N = length(X); % number of training examples
	S = length(hmm.t); % number of states

	Nt = zeros(S, 1);
	NT = zeros(S, S);
	means = zeros(S, 1);
	vars = zeros(S, 1);
	W = zeros(S, 

	%% E-Step
	for i=1:N
		x = X{i};	% particular vector of emissions (DATA)
		m = length(x); % number of timesteps in hidden chain
		log_a = forward(x, hmm);
		log_b = backward(x, hmm);
		loglik = loglik + log_sum_exp(log_a(:,end),1); % add log Pr(DATA)

		% update initial state counts	
    Nt = Nt + norm_exp(log_a(:,1)+log_b(:,1));

		% update transition counts
    for t=1:m-1	
			log_A = repmat(log_a(:,t),1,k);
      log_B = repmat(log_b(:,t+1),1,k);
			log_E = repmat(hmm.logEall, 1, k);
			NT = NT + norm_exp(log_A+hmm.log_T+log_E'+log_B'); 
    end

		
	end % end LOOP over training examples

	%% M-Step
	hmm.t = Nt / sum(Nt);	
	hmm.log_t = log(hmm.t);

	hmm.T = NT ./ repmat(sum(NT, 2), 1, S);
	hmm.log_T = log(hmm.T);
end

%function [hmm, loglik] = one_example(x, hmm)

%end
