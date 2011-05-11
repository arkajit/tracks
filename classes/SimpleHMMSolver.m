classdef SimpleHMMSolver < handle

	properties
		Dmax
		Vmax
		nBins
		odds
		hmm
	end

	methods

		function [self] = SimpleHMMSolver(Dmax, Vmax, nBins, odds)
			self.Dmax = Dmax;
			self.Vmax = Vmax;

			if (nargin < 3)
				self.nBins = 20;
			else
				self.nBins = nBins;
			end

			if (nargin < 4)
				self.odds = 10;
			else 
				self.odds = odds;
			end
		end

		function [self] = selectModel(self, tau)
			maxstd = sqrt(2*self.Dmax*tau);
			maxmean = self.Vmax*tau;

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

			S = self.nBins^2;
			mc = MarkovChain.fromOdds(S, self.odds);
			means = reshape(repmat(mus, self.nBins, 1), S, 1);
			stddevs = reshape(repmat(sigmas', self.nBins, 1), S, 1);
			self.hmm = CHMM(mc, means, stddevs);
		end

		function [D, V] = infer(self, track)
			D = zeros(size(track.V)); % 3D
			V = zeros(size(track.V));

			for i=1:3
				states = self.hmm.viterbi(track.steps(:,i));
				V(:,i) = self.hmm.means(states) ./ track.tau;
				D(:,i) = (self.hmm.stddevs(states) .^ 2) ./ (2 * track.tau);
			end
		end

		function [D, V] = solve(self, track)
			self.selectModel(track.tau);
			[D, V] = self.infer(track);
		end

		function [D, V, errs] = test(self, tracks)
			N = length(tracks);
			if (~N)
				return;
			else
				tau = tracks(1).tau;
				self.selectModel(tau);
			end
			D = cell(N, 1);
			V = cell(N, 1);
			errs = zeros(N, 4);

			for i=1:N
				t = tracks(i);
				[D{i}, V{i}] = self.infer(t);
				errs(i,:) = t.compare(D{i}(:,1), V{i});
			end
		end

	end

end
