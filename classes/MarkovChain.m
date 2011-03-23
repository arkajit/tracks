classdef MarkovChain
    %MARKOVCHAIN A markov chain
    
    properties
      S         % number of states
      start     % initial distribution of states
      trans     % S x S transition matrix
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

		end
    
end
