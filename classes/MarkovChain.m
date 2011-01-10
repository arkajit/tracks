classdef MarkovChain
    %MARKOVCHAIN A markov chain
    %   Detailed explanation goes here
    
    properties
      K         % number of states
      start     % initial distribution of states
      trans     % K x K transition matrix
    end
    
    methods
    
      function self = MarkovChain(t, T)
        self.K = size(t,1);
        if (size(T) ~= [self.K, self.K])
          disp(Error: transition matrix has wrong dimensions)
          return;
        end
        self.start = t;
        self.trans = T;
      end
    end
    
end
