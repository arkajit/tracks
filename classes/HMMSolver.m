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

		% Initialize a HMM with random parameters
		% Will be used to start EM
		%
		%	@param 	tau		double	timestep in seconds
		% @return hmm		CHMM		a random, initial model
		function [hmm] = init(self, tau)
      maxstd = sqrt(2*self.Dmax*tau);
      maxmean = self.Vmax*tau;
			hmm = CHMM.random(self.S, maxmean, maxstd);
		end

	end

	methods (Static)

		function [errs] = errors(n, A, B)
			errs = -1;
			for f=perms(1:n)
				nErrs = 0;
				
				% compute the errors with this permutation
				for i=1:length(A)
					nErrs = nErrs + sum(A{i} ~= f(B{i}));
					if (errs >= 0 && nErrs > errs)
						break;
					end
				end

				% update best error count so far
				if (errs < 0 || nErrs < errs)
					errs = nErrs;
				end
			end
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
		% @param maxIter	int					number of EM iterations
		function [self] = train(self, tracks, maxIter)
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
			self.hmm = self.hmm.em(X, maxIter);
		end

		function [D, V] = test(self, tracks)
			N = length(tracks);		% no. of test examples
			D = cell(N, 1);
			V = cell(N, 1);

			for i=1:N
				[D{i}, V{i}] = self.infer(tracks{i});
			end
		end

	end
end
