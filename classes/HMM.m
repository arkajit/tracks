classdef HMM

	properties (Abstract)
		S				% number of states
		pi			% initial state distribution
		M				% state transition matrix
	end

	properties (Dependent)
		log_pi
		log_M
	end
	
	methods (Abstract)
		[prob] = E(self, s, x);		% emission function
		[x] = emit(self, s);			% sample the emission distribution
	end

	methods 
		function [lpi] = get.log_pi(self)
			lpi = log(self.pi); 
		end

		function [lm] = get.log_M(self)
			lm = log(self.M);
		end
	
		function [prob] = log_E(self, s, x)
			prob = log(self.E(s, x));
		end
		
		function [probs] = Es(self, x)
			probs = zeros(self.S, 1);
			for s=1:self.S
				probs(s) = self.E(s, x);
			end	
		end

		function [probs] = log_Es(self, x)
			probs = log(self.Es(x));
		end
	end

	methods

		% sample N sequences of observations of length T from this HMM
		function [X, S] = sample(self, N, T)
			X = zeros(T, N);
			S = zeros(T, N);
			for i=1:N
				S(1,i) = mnsmpl(self.pi);
				X(1,i) = self.emit(S(1,i));
				
				for t=2:T
					S(t,i) = mnsmpl(self.M(S(t-1,i),:));
					X(t,i) = self.emit(S(t,i));
				end
			end

			function [r] = mnsmpl(p) 
				p = reshape(p, 1, length(p));
				r = find(mnrnd(1, p)==1);
			end
		end

		function [S] = infer(self, X)
			N = size(X,2);
			S = zeros(size(X));
			for i=1:N
				S(:,i) = self.viterbi(X(:,i));
			end
		end

		function [state, log_v] = viterbi(self, x) 
			S = self.S;						 % number of states
			T = length(x);         % length of sequence
			log_v = zeros(S, T); 

			for s=1:S 	% for each state
				log_v(s,1) = self.log_pi(s) + self.log_E(s,x(1)); 
			end

			for t=2:T
				log_V = repmat(log_v(:,t-1),1,S); % each col is log_v(:,t-1)

				best = max(log_V+self.log_M,[],1)';
				for s=1:S 
					log_v(s,t) = best(s) + self.log_E(s,x(t));
				end
			end

			% backtrack
			state = zeros(T, 1);
			[tmp,ind] = max(log_v(:,T));
			state(T) = ind(1); 
			for t=T-1:-1:1,
					[tmp,ind] = max(log_v(:,t)+self.log_M(:,state(t+1))); 
					state(t) = ind(1); 
			end
		end

		function [log_a] = forward(self, x)
			S = self.S;				% number of states
			T = length(x);		% length of sequence
			log_a = zeros(S, T);
		
			for s=1:S
				log_a(s, 1) = self.log_pi(s) + self.log_E(s,x(1));
			end

			for t=2:T
				log_A = repmat(log_a(:,t-1),1,S);	% each col is log_a(:,t-1)
				prev = log_sum_exp(log_A+self.log_M, 1)';

				for s=1:S
					log_a(s,t) = prev(s) + self.log_E(s,x(t));
				end
			end
		end

		function [log_b] = backward(self, x)
			S = self.S;				% number of states
			T = length(x);		% length of sequence
			log_b = zeros(S, T);
			% Init: b(:,T) = 1 => log_b(:,T) = 0

			for t=T-1:-1:1
				log_e = zeros(S, 1);
				for s=1:S
					log_e(s) = self.log_E(s, x(t+1));
				end

				log_E = repmat(log_e, 1, S);					% SxS
				log_B = repmat(log_b(:,t+1), 1, S);		% SxS
				log_b(:,t) = log_sum_exp(self.log_M + log_E' + log_B', 2);
			end
		end

	end

	methods (Static)

		function [log_p] = log_sum_exp(log_P,ind)
			log_max  = max(log_P,[],ind); 
			log_Pmax = repmat(log_max,size(log_P)./size(log_max));
			I = find(log_Pmax>-Inf); % shouldn't subtract -Inf 
			log_Pmax(I) = log_P(I) - log_Pmax(I); 
			log_p = log_max + log(sum(exp(log_Pmax),ind));
		end

		function [p] = norm_exp(log_p) 
			p = exp(log_p(:)-max(log_p(:))); 
			p = reshape(p./sum(p), size(log_p) );
		end

	end
end
