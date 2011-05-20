classdef CHMM < HMM
	%CHMM a continuous hidden markov model

	properties
		chain			% MarkovChain
		means 		% means of states
		stddevs 	% standard deviations of states	
		options		% learning options
		S
	end

	properties (Dependent)
		pi
		M
	end

	methods 			% property getters
		%function S = get.S(self)
	%		S = self.chain.S;
	%	end

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
			if (nargin < 2 || ~isnumeric(odds))
				odds = 10;
			end
			hmm.chain = MarkovChain.fromOdds(S, odds);
		end

		function [hmm, L, hmms] = fit(S, X, options, nRestarts, maxIter)
			if (nargin < 3)
				options = [];
			end
	
			if (nargin < 4)
				nRestarts = 2;
			end

			if (nargin < 5)
				maxIter = 15;	
			end

			for i=1:nRestarts
				fprintf('EM restart %d: \n', i);
				if (isstruct(options))
					hmm0 = CHMM.fromOdds(S);
					hmm0.options = options;
				else
					hmm0 = CHMM.random(S);
				end
				[hmm1, L1] = hmm0.em(X, maxIter);	
				hmms(i) = hmm1;
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
			self.S = mc.S;
			[self.means, order] = sort(means);
			self.stddevs = stddevs(order);	
			if (length(self.means) ~= self.chain.S || ...
					(length(self.stddevs) ~= self.chain.S))
				disp('Error: incorrect number of means or stddevs');
			end

			if (nargin < 4)
				self.options.learnStart = true;
				self.options.learnTrans = true;
				self.options.learnOdds = false;
				self.options.learnEmit = true;
			else
				self.options = options;
			end
		end

		% emission matrix: assumes x is a column vector
		function probs = Es(self, x)
			T = length(x);
			probs = zeros(self.S, T);
			probs = normpdf(repmat(x', self.S, 1), ...
											repmat(self.means, 1, T), ...
											repmat(self.stddevs, 1, T));
			probs(isnan(probs)) = 0;	% non-positive stddevs will give NaNs
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
		% @param 		maxIter	int					max number of EM iterations (default: 15)
		% @return 	hmm 		CHMM
		% @return 	L				double			log-likelihood
		function [hmm, L, logliks] = em(self, X, maxIter)
			if (nargin < 3)
				maxIter = 15;
			end
			logliks = nan(maxIter, 1);		

			fprintf('Starting EM...\n');
			[hmm, L] = self.one_em_iter(X);
			logliks(1) = L;
			fprintf('Initial loglik = %f\n', L);

			L0 = L-abs(L);
			iter = 1;
            % immediately stop if log likelihood ever decreases
			while (iter < maxIter && L-L0 > abs(L)*1e-6)
				L0 = L;
				hmm0 = hmm;
				[hmm, L] = hmm.one_em_iter(X);
				if (L-L0 < 0)
					fprintf('**');
				end
				fprintf('Iter %d loglik = %f\n', iter, L);
				iter = iter + 1;
				logliks(iter) = L;
			end
			
			if (L-L0 < 0)
				L = L0;
				hmm = hmm0;
				fprintf('Loglik decreased. Revert to last loglik.\n');
			end
			fprintf('Final loglik = %f\n', L);
		end

	end

	methods (Access=private)

		% Perform one iteration of EM starting from current model (self).
		% Return a new model with increased likelihood. If there are too many
		% states in the model, may eliminate extra states by having NaNs for means
		% and stddevs.
		%
		% @param 	X 			mat		T x N		N training examples of length T
		%	@return hmm			CHMM
		% @return loglik 	double						the log likelihood of the new model
		function [hmm, loglik] = one_em_iter(self, X)
			S = self.S; 						% number of states
			N = size(X, 2);					% number of training examples
			T = size(X, 1);					% number of timesteps

			% Laplace correction: start with a pseudo-count of 1 instead of 0
			Npi = ones(S, 1);			% counts of initial states
			NM = ones(S, S);				% counts of transitions

			NW = cell(N, 1);				% unnormalized weights of each state
			PS = zeros(S, 1);				% probs of states across all data
			PM = zeros(S, 1);				% unnormalized counts of means
			PV = zeros(S, 1);				% unnormalized counts of variances

			means = zeros(S, 1);		% estimates of means
			stddevs = zeros(S, 1);	% estimates of standard deviations
			loglik = 0;							% the log likelihood of the new model

			%% E-Step
			for i=1:N
				%% compute initial statistics
				x = X(:,i);
				log_a = self.forward(x);
				log_b = self.backward(x);
				probX = logsumexp(log_a(:,end),1);
				if (isnan(probX))
					fprintf('Error: probabilities are not numbers!\n');
				end
				loglik = loglik + probX;

				%% calculate update weights
				if (self.options.learnEmit)
					log_w = log_a + log_b - probX;
					W = exp(log_w);
					NW{i} = W;

					PS = PS + sum(W, 2);
					PM = PM + W * x;
				end

				%% update initial state counts	
				if (self.options.learnStart) 
					Npi = Npi + norm_exp(log_a(:,1)+log_b(:,1));
				end

				%% update transition counts
				if (self.options.learnTrans)
					log_e = log(self.Es(x));
					log_M = log(self.M);
					for t=1:T-1	
						log_A = repmat(log_a(:,t),1,S);
						log_B = repmat(log_b(:,t+1),1,S);
						log_E = repmat(log_e(:,t+1),1,S);
						NM = NM + norm_exp(log_A+log_M+log_E'+log_B'); 
					end
				end

			end	% END loop over training examples

			%% Compute means and variance
			if (self.options.learnEmit)
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
			end

			%% M-Step (for pi and M)
			if (self.options.learnStart)
				pi = Npi / sum(Npi);
			else
				pi = self.chain.start;
			end

			odds = NaN;
			if (self.options.learnTrans)
				M = NM ./ repmat(sum(NM, 2), 1, S);
			elseif (self.options.learnOdds)
				stay = diag(NM);
				swap = NM - diag(stay);
				odds = sum(stay) / sum(sum(swap));
				M = MarkovChain.oddsMatrix(S, odds);
			else
				M = self.chain.trans;
			end

			hmm = CHMM(MarkovChain(pi, M), means, stddevs, self.options);
			hmm.chain.odds = odds;
		end

	end

end
