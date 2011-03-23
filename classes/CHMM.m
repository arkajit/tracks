classdef CHMM < HMM
	%CHMM a continuous hidden markov model

	properties
		chain			% MarkovChain
		means 		% means of states
		stddevs 	% standard deviations of states	
	end

	properties (Dependent)
		S
		pi
		M
	end

	methods 			% property getters
		function S = get.S(self)
			S = self.chain.S;
		end

		function pi = get.pi(self)
			pi = self.chain.start;
		end

		function M = get.M(self)
			M = self.chain.trans;
		end
	end

	methods

		function self = CHMM(mc, means, stddevs)
			self.chain = mc;
			self.means = means;
			self.stddevs = stddevs;	
			if (length(self.means) ~= self.chain.S || ...
					(length(self.stddevs) ~= self.chain.S))
				disp('Error: incorrect number of means or stddevs');
			end
		end

		% emission function
		function prob = E(self, s, x)
			if (self.stddevs(s))
				prob = normpdf(x, self.means(s), self.stddevs(s));
			else
				prob = double(x == means(s));
			end	
		end

		function hmm = em(self, x)
			[hmm, L] = self.one_em_iter(x);
			disp(sprintf('Initial loglik = %f', L));

			L0 = L-abs(L);
			while (L-L0 > abs(L)*1e-6)
				L0 = L;
				[hmm, L] = hmm.one_em_iter(x);
			end

			disp(sprintf('Final loglik = %f', L));
		end

	end

	methods (Access=private)

		function [hmm, loglik] = one_em_iter(self, x)
			S = self.S; 					% number of states
			T = length(x); 				% number of timesteps in hidden chain

			Npi = zeros(S, 1);			% counts of initial states
			NM = zeros(S, S);			% counts of transitions
			W = zeros(S, T);			% weights to update mean and variances

			%% E-Step
			log_a = self.forward(x);
			log_b = self.backward(x);
			loglik = log_sum_exp(log_a(:,end),1); 		% log Pr(DATA)

			% calculate update weights
			log_w = log_a + log_b - loglik;
			W = exp(log_w);
			W = W ./ repmat(sum(W, 2), 1, T); 				% normalize weights

			% update initial state counts	
			Npi = Npi + norm_exp(log_a(:,1)+log_b(:,1));

			% update transition counts
			for t=1:T-1	
				log_A = repmat(log_a(:,t),1,S);
				log_B = repmat(log_b(:,t+1),1,S);
				log_E = repmat(self.log_Es(x(t+1)),1,S);
				NM = NM + norm_exp(log_A+self.log_M+log_E'+log_B'); 
			end

			%% M-Step
			pi = Npi / sum(Npi);
			M = NM ./ repmat(sum(NM, 2), 1, S);
			means = W * x;
			devs = repmat(x, 1, S) - repmat(means', T, 1);		% TxS: x_t - mu_s
			stddevs = sqrt(diag(W * (devs.^2)));
			hmm = CHMM(MarkovChain(pi, M), means, stddevs);
		end	

	end

end
