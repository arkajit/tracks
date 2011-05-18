classdef HMM
	properties (Abstract)
		S				% number of states
		pi			% initial state distribution
		M				% state transition matrix
	end

	methods (Abstract)
		[E] = Es(self, x);				% emission matrix (S x T) where T = length(x)
		[x] = emit(self, s);			% sample the emission distribution
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
			log_e = log(self.Es(x));
			log_M = log(self.M);

			log_v = zeros(S, T); 
			log_v(:,1) = log(self.pi) + log_e(:,1);	

			for t=2:T
				%log_V = repmat(log_v(:,t-1),1,S); % each col is log_v(:,t-1)
				best = max(bsxfun(@plus, log_v(:,t-1), log_M),[],1)';
				log_v(:,t) = best + log_e(:,t);
			end

			% backtrack
			state = zeros(T, 1);
			[tmp,ind] = max(log_v(:,T));
			state(T) = ind(1); 
			for t=T-1:-1:1,
					[tmp,ind] = max(log_v(:,t)+log_M(:,state(t+1))); 
					state(t) = ind(1); 
			end
		end

		function [log_a] = forward(self, x)
			S = self.S;				% number of states
			T = length(x);		% length of sequence
			log_e = log(self.Es(x));
			log_M = log(self.M);

			log_a = zeros(S, T);
			log_a(:,1) = log(self.pi) + log_e(:,1);

			for t=2:T
				%log_A = repmat(log_a(:,t-1),1,S);	% each col is log_a(:,t-1)
				prev = logsumexp(bsxfun(@plus, log_a(:,t-1), log_M), 1)';
				log_a(:,t) = prev + log_e(:,t);
			end
		end

		function [log_b] = backward(self, x)
			S = self.S;				% number of states
			T = length(x);		% length of sequence
			log_e = log(self.Es(x));
			log_M = log(self.M);

			log_b = zeros(S, T);
			% Init: b(:,T) = 1 => log_b(:,T) = 0

			for t=T-1:-1:1
				%log_E = repmat(log_e(:,t+1), 1, S);		% SxS
				%log_B = repmat(log_b(:,t+1), 1, S);		% SxS
				log_a = log_e(:,t+1) + log_b(:,t+1);
				%log_b(:,t) = logsumexp(log_M + log_E' + log_B', 2);
				log_b(:,t) = logsumexp(bsxfun(@plus, log_a', log_M), 2);
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

		% assumes we're comparing vector of discrete state labels
		function [frac, errs, order] = errors(n, A, B)
			errs = -1;
			[T, alen] = size(A);
			order = [];
			for f=perms(1:n)'
				nErrs = 0;
				
				% compute the errors with this permutation
				for i=1:alen
					nErrs = nErrs + sum(A(:,i) ~= f(B(:,i)));
					if (errs >= 0 && nErrs > errs)
						break;
					end
				end

				% update best error count so far
				if (errs < 0 || nErrs < errs)
					errs = nErrs;
					order = f;
				end
			end

			frac = errs / (alen * T);
		end

	end
end
