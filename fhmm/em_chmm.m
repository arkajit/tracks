function [hmm] = em_chmm(x, hmm)
% em_chmm - Expectation Maximization for a continuous Gaussian HMM
%		x - mx1 vector of observations (a single training example)
%		hmm - initial guess for the parameters
%			requires t, T, means, stddevs fields to be set to init values
%			returns the MLE of hmm parameters including t, T, means, stddevs

[hmm, L] = one_em_iter(x, hmm);
disp(sprintf('Initial loglik = %f', L));

L0 = L-abs(L);
while (L-L0 > abs(L)*1e-6)
	L0 = L;
	[hmm, L] = one_em_iter(x, hmm);
end

hmm.loglik = L;
disp(sprintf('Final loglik = %f', L));

% -----------------------------------------
function [hmm, loglik] = one_em_iter(x, hmm)
	S = length(hmm.t); 		% number of states
	m = length(x); 				% number of timesteps in hidden chain

	Nt = zeros(S, 1);			% counts of initial states
	NT = zeros(S, S);			% counts of transitions
	W = zeros(S, m);			% weights to update mean and variances

	%% E-Step
	log_a = forward(x, hmm);
	log_b = backward(x, hmm);
	loglik = log_sum_exp(log_a(:,end),1); 		% log Pr(DATA)

	% calculate update weights
	log_w = log_a + log_b - loglik;
	W = exp(log_w);
	W = W ./ repmat(sum(W, 2), 1, m); 				% normalize weights

	% update initial state counts	
	Nt = Nt + norm_exp(log_a(:,1)+log_b(:,1));

	% update transition counts
	for t=1:m-1	
		log_A = repmat(log_a(:,t),1,k);
		log_B = repmat(log_b(:,t+1),1,k);
		log_E = repmat(hmm.logEall, 1, k);
		NT = NT + norm_exp(log_A+hmm.log_T+log_E'+log_B'); 
	end

	%% M-Step
	hmm.t = Nt / sum(Nt);	
	hmm.log_t = log(hmm.t);

	hmm.T = NT ./ repmat(sum(NT, 2), 1, S);
	hmm.log_T = log(hmm.T);

	hmm.means = W * x;
	hmm.stddevs = sqrt(W * ((x - hmm.means).^2));
