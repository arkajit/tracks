classdef HMMSolver < TrackSolver
	% represents the logic used to solve a Track using an HMM

	properties
		Dmax
		Vmax
		S 		% expected number of states
		hmm
	end	

	methods (Access = private)
	
		% Using selected HMM, perform inference to find the most likely sequence of
		% motion parameters that could generate this track. Computes a Viterbi
		% decoding.
		% 
		% @param track Track
		% 
		% @return D mat
		% @return V mat 
		function [D, V] = infer(self, track)
      D = zeros(size(track.V)); % full 3D (should be similar in each dimension)
      V = zeros(size(track.V));

      for i=1:3
				states = self.hmm.viterbi(track.steps(:,i));
				V(:,i) = self.hmm.means(states) ./ track.tau;
				D(:,i) = (self.hmm.stddevs(states) .^ 2) ./ (2 * track.tau);
      end
		end	

		function [err] = computeError(self, errD, errV) 
      errD_scaled = errD / self.Dmax;
      errV_scaled = errV ./ (2*self.Vmax);
			
			% max error norm is 2
      % since each of the four components has max error 1
      err = norm([errV_scaled; errD_scaled]) / 2;  
		end

		% Initialize a HMM with random parameters
		% Will be used to start EM
		%
		%	@param 	tau		double	timestep in seconds
		% @return hmm		CHMM		a random, initial model
		function [hmm] = init(self, tau)
      maxstd = sqrt(2*self.Dmax*tau);
      maxmean = self.Vmax*tau;
			sigs = rand(self.S, 1) * maxstd;	
			mus = rand(self.S, 1) .* sign(randn(self.S, 1)) * maxmean;
			mc = MarkovChain.random(self.S);
			hmm = CHMM(mc, mus, sigs);
		end

	end

	methods
		
		% constructor
		function [self] = HMMSolver(Dmax, Vmax, S) 
			self.Dmax = Dmax;
			self.Vmax = Vmax;
			self.S = S;
		end

		% Select the ML HMM that explains the observed tracks.
		% Require all tracks to have the same timestep.
		% 
		% @param tracks 	cellarray 	Nx1
		function [self] = train(self, tracks)
			N = length(tracks);
			if (~N)
				return;
			end

			X = cell(3*N, 1);
			tau = 0;
			for i=1:N
				track = tracks{i};

				if (~tau) 
					tau = track.tau;
				elseif (tau ~= track.tau)
					disp('Error: tracks must have the same timestep');
					return;
				end

				for j=1:3
					X{3*(i-1)+j} = track.steps(:,j);
				end
			end

			self.hmm = self.init(tau);
			self.hmm = self.hmm.em(X);
		end

		% Solve a track for its underlying motion parameters using a factorial HMM.
		% 
		% @param track Track
		function [D, V, err] = solve(self, track)
			[D, V] = self.infer(track);
			[errD, errV] = track.compare(D, V);
			err = self.computeError(errD, errV);	
		end

		function [D, V] = solveAll(self, tracks)
			N = length(tracks);		% no. of test examples
			D = cell(N, 1);
			V = cell(N, 1);

			for i=1:N
				[D{i}, V{i}] = self.infer(tracks{i});
			end
		end

	end
end
