classdef HMMSolver < TrackSolver
	% represents the logic used to solve a Track using an HMM

	properties
		Dmax
		Vmax
		nBins
		odds
		hmm
		hmms
	end	

	methods (Access = private)
	
		function [self] = selectModel(self, track) 
      maxstd = sqrt(2*self.Dmax*track.tau);
      maxmean = self.Vmax*track.tau;
    
      if self.Dmax == 0
        sigmas = [0.1]; % disallowing all noise is too error-prone
      else 
        sigincr = maxstd/self.nBins;
        sigbins = 0:sigincr:(maxstd-sigincr);
				
				% just get the centers by shifting halfbin
        sigmas = sigbins + sigincr/2; 
      end

      if self.Vmax == 0
        mus = [0];
      else
        muincr = (2*maxmean)/self.nBins;
        mubins = -maxmean:muincr:(maxmean-muincr);

				% now we can never select V=0
        mus = mubins + muincr/2; 
      end

			self.hmm = self.initFHMM(mus, sigmas);

			% learn three different HMMs for each direction
			% using the same initial guess supplied by the odds
			self.hmms = CHMM.empty(3, 0);
			for i=1:3
				%self.hmms(i) = self.hmm.em(track.steps(:,i));
				self.hmms(i) = self.hmm; % Temporary since EM still buggy. 
			end
		end

		% Using selected HMM, perform inference to find the most likely sequence of
		% motion parameters that could generate this track. Computes a Viterbi
		% decoding.
		% 
		% @param track Track
		% 
		% @return D vec
		% @return V mat 
		function [D, V] = infer(self, track)
      D = zeros(size(track.D));
      D3 = zeros(size(track.V)); % want 3D
      V = zeros(size(track.V));

      for i=1:3
				hmm = self.hmms(i);
				states = hmm.viterbi(track.steps(:,i));

				V(:,i) = hmm.means(states) ./ track.tau;
				D3(:,i) = (hmm.stddevs(states) .^ 2) ./ (2 * track.tau);
      end
		
			% get average of x,y,z predictions
      D = mean(D3, 2);
		end	

		function [err] = computeError(self, errD, errV) 
      errD_scaled = errD / self.Dmax;
      errV_scaled = errV ./ (2*self.Vmax);
			
			% max error norm is 2
      % since each of the four components has max error 1
      err = norm([errV_scaled; errD_scaled]) / 2;  
		end


		% Initialize a fHMM from provided odds parameters
		% Will be used as starting point for EM
		function [hmm] = initFHMM(self, mus, sigs)
			P = length(mus);
			Q = length(sigs);
			S = P*Q;
			pi = ones(S, 1);
			if S == 1
				M = [1];
			else
				M = ones(S) + ((S-1)*self.odds - 1)*eye(S); % diag. entries = odds*(S-1)
			end
			
			mc = MarkovChain(pi, M);
			means = reshape(repmat(mus, Q, 1), S, 1);
			stddevs = reshape(repmat(sigs', P, 1), S, 1);
			hmm = CHMM(mc, means, stddevs);
		end

	end

	methods
		
		% constructor
		function [self] = HMMSolver(Dmax, Vmax, nBins, odds) 
			self.Dmax = Dmax;
			self.Vmax = Vmax;
			self.nBins = nBins;
			if (nargin == 4)
				self.odds = odds;
			else
				self.odds = 10;
			end
		end

		% Solve a track for its underlying motion parameters using a factorial HMM.
		% 
		% @param track Track
		function [D, V, err] = solve(self, track)
			self = self.selectModel(track);
			[D, V] = self.infer(track);
			[errD, errV] = track.compare(D, V);
			err = self.computeError(errD, errV);	
		end

	end
end
