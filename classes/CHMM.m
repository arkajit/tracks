classdef CHMM < HMM
	%CHMM a continuous hidden markov model

	properties
		chain			% MarkovChain
		means 		% means of states
		stddevs 	% standard deviations of states	
		options		% learning options
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

	methods (Static)
	
		function hmm = random(S, mumax, stdmax)
			if (nargin == 1)
				mumax = 10;
				stdmax = 2;
			end
			mc = MarkovChain.random(S);
			means = rand(S, 1) .* sign(randn(S, 1)) * mumax;
			stddevs = rand(S, 1) * stdmax;
			hmm = CHMM(mc, means, stddevs);
		end

		function hmm = fromOdds(S, odds)
			hmm = CHMM.random(S);
			hmm.chain = MarkovChain.fromOdds(S, odds);
		end

		function [hmm, L] = fit(S, X, nRestarts, maxIter)
			if (nargin < 3)
				nRestarts = 2;
			end

			if (nargin < 4)
				maxIter = 10;	
			end

			for i=1:nRestarts
				disp(sprintf('EM restart %d: ', i));
				hmm0 = CHMM.random(S);
				[hmm1, L1] = hmm0.em(X, maxIter);	
				if (i == 1 || isnan(L) || L1 > L)
					hmm = hmm1;
					L = L1;
				end
			end
		end

	end

	methods

		function self = CHMM(mc, means, stddevs, options)
			self.chain = mc;
			self.means = means;
			self.stddevs = stddevs;	
			if (length(self.means) ~= self.chain.S || ...
					(length(self.stddevs) ~= self.chain.S))
				disp('Error: incorrect number of means or stddevs');
			end

			if (nargin < 4)
				self.options.learnStart = true;
				self.options.learnTrans = true;
			else
				self.options = options;
			end
		end

		% emission function
		function prob = E(self, s, x)
			prob = normpdf(x, self.means(s), self.stddevs(s));
			if isnan(prob)
				prob = double(x == self.means(s));
			end	
		end

		% sampling the emission function
		function x = emit(self, s)
			x = self.means(s) + randn * self.stddevs(s);
		end

		% Overlay and plot all the normal emission distributions for given range.
		%
		% @param 	x 	vec
		function plotNormals(self, x)
			if (nargin == 1)
				x = -15:0.01:15;
			end
			y = zeros(length(x), self.S);
			for i=1:self.S
				y(:, i) = normpdf(x, self.means(i), self.stddevs(i));
			end
			plot(x, y);
		end

		% Perform Expectation Maximization on a set of training examples.
		%	Returns the HMM that maximizes the likelihood of the observed data.
		% 
		%	@param 		X 			mat					T x N
		% @param 		maxIter	int					max number of EM iterations (default: 10)
		% @return 	hmm 		CHMM
		% @return 	L				double			log-likelihood
		function [hmm, L] = em(self, X, maxIter)
			if (nargin < 3)
				maxIter = 10;
			end

			disp('Starting EM...');
			[hmm, L] = self.one_em_iter(X);
			disp(sprintf('Initial loglik = %f', L));

			L0 = L-abs(L);
			iter = 1;
			while (iter < maxIter && L-L0 > abs(L)*1e-6)
				L0 = L;
				[hmm, L] = hmm.one_em_iter(X);
				disp(sprintf('Iter %d loglik = %f', iter, L));
				iter = iter + 1;
			end

			disp(sprintf('Final loglik = %f', L));
		end

	end

	methods (Access=private)

		% Perform one iteration of EM starting from current model (self).
		% Return a new model with increased likelihood.
		%
		% @param 	X 			mat		T x N		N training examples of length T
		%	@return hmm			CHMM
		% @return loglik 	double						the log likelihood of the new model
		function [hmm, loglik] = one_em_iter(self, X)
			S = self.S; 						% number of states
			N = size(X, 2);					% number of training examples
			T = size(X, 1);					% number of timesteps

			Npi = zeros(S, 1);			% counts of initial states
			NM = zeros(S, S);				% counts of transitions

			NW = cell(N, 1);				% unnormalized weights of each state
			PS = zeros(S, 1);				% probs of states across all data
			PM = zeros(S, 1);				% unnormalized counts of means
			PV = zeros(S, 1);				% unnormalized counts of variances

			means = zeros(S, 1);		% estimates of means
			stddevs = zeros(S, 1);	% estimates of standard deviations
			loglik = 0;							% the log likelihood of the new model

			%% E-Step
			for i=1:N
				x = X(:,i);
				W = zeros(S, T);			% weights to update mean and variances

				log_a = self.forward(x);
				log_b = self.backward(x);
				probX = log_sum_exp(log_a(:,end),1);
				loglik = loglik + probX;

				% calculate update weights
				log_w = log_a + log_b - probX;
				W = exp(log_w);
				NW{i} = W;

				PS = PS + sum(W, 2);
				PM = PM + W * x;

				if (self.options.learnStart) 
					% update initial state counts	
					Npi = Npi + norm_exp(log_a(:,1)+log_b(:,1));
				end

				if (self.options.learnTrans)
					% update transition counts
					for t=1:T-1	
						log_A = repmat(log_a(:,t),1,S);
						log_B = repmat(log_b(:,t+1),1,S);
						log_E = repmat(self.log_Es(x(t+1)),1,S);
						NM = NM + norm_exp(log_A+self.log_M+log_E'+log_B'); 
					end
				end

			end	% END loop over training examples

			%% Compute means and variance
			means = PM ./ PS;
			for i=1:N
				x = X(:,i);
				T = length(x);
				if (~T)
					continue;
				end

				W = NW{i};
				for s=1:S
					PV(s) = PV(s) + W(s,:) * (x-means(s)).^2;
				end
			end
			stddevs = sqrt(PV ./ PS);

			%% M-Step (for pi and M)
			if (self.options.learnStart)
				pi = Npi / sum(Npi);
			else
				pi = self.chain.start;
			end

			if (self.options.learnTrans)
				M = NM ./ repmat(sum(NM, 2), 1, S);
			else
				M = self.chain.trans;
			end

			hmm = CHMM(MarkovChain(pi, M), means, stddevs, self.options);
		end

	end

end
