classdef MarkovChain
    %MARKOVCHAIN A markov chain
    
    properties
      S         % number of states
      start     % initial distribution of states
      trans     % S x S transition matrix
			odds			% if the transition matrix is parameterized by odds
    end
    
    methods
    
      function self = MarkovChain(pi, M)
        self.S = size(pi,1);
        if (size(M) ~= [self.S, self.S])
          disp('Error: transition matrix has wrong dimensions')
          return;
        end
        self.start = pi ./ sum(pi);
        self.trans = M ./ repmat(sum(M, 2), 1, self.S);
				self.odds = NaN;
      end

			function [states] = sample(self, T)
				states = zeros(T, 1);
				states(1) = mnsmpl(self.start);
				for t=2:T
					states(t) = mnsmpl(self.trans(states(t-1),:));
				end

				function [r] = mnsmpl(p) 
					p = reshape(p, 1, length(p));
					r = find(mnrnd(1, p)==1);
				end
			end

    end

		methods (Static)
			
			function mc = uniform(S)
				pi = ones(S, 1);
				M = ones(S, S);
				mc = MarkovChain(pi, M);
			end

			function mc = random(S)
				pi = rand(S, 1);
				M = rand(S, S);
				mc = MarkovChain(pi, M);
			end

			function mc = fromOdds(S, odds)
				pi = rand(S, 1);
				M = MarkovChain.oddsMatrix(S, odds);
				mc = MarkovChain(pi, M);
			end

			function M = oddsMatrix(S, odds)
				M = ones(S, S);
				if (S > 1)
					for i=1:S
						M(i,i) = (S-1)*odds; 
					end
				end
			end

		end
    
end
