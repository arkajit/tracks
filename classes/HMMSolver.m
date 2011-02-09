classdef HMMSolver
	% represents the logic used to solve a Track using an HMM

	properties
		Dmax
		Vmax
		nBins
		hmm
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

			% uniform prior, 10:1 odds against switching	
			self.hmm = fact_hmm(mus, sigmas, [], [], 10, 10);
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
        states = viterbi(track.steps(:,i), self.hmm);
        [Si, inds] = self.hmm.to_pairs(states);

        V(:,i) = Si(:,1) ./ track.tau;
        D3(:,i) = (Si(:,2) .^ 2) ./ (2 * track.tau);
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

	end

	methods
		
		% constructor
		function [self] = HMMSolver(Dmax, Vmax, nBins) 
			self.Dmax = Dmax;
			self.Vmax = Vmax;
			self.nBins = nBins;
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
