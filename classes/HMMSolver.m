classdef HMMSolver < TrackSolver
	% represents the logic used to solve a Track using an HMM

	properties
		S 				% expected number of states
		hmms			% list of hmms for each dimension
		options		% learning options (default, if not provided)
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
				if (~ isempty(self.hmms))
					hmm = self.hmms(i);
				else 
					hmm = self.hmm;
				end

				states = hmm.viterbi(track.steps(:,i));
				V(:,i) = hmm.means(states) ./ track.tau;
				D(:,i) = (hmm.stddevs(states) .^ 2) ./ (2 * track.tau);
      end
		end	

		% Initialize a HMM with random parameters
		% Will be used to start EM
		%
		%	@param 	tau		double	timestep in seconds
		% @return hmm		CHMM		a random, initial model
		function [hmm] = init(self, tau)
			hmm = CHMM.random(self.S);
		end

	end

	methods (Static)

		function [frac] = pctErr(A, B)
			frac = 100 * norm(A-B) / norm(A);
		end

		function [X, Y, Z, tau] = readTracks(tracks)
			N = length(tracks);
			if (~N)
				return;
			end
			T = tracks(1).T - 1;

			X = zeros(T,N);
			Y = zeros(T,N);
			Z = zeros(T,N);

			for i=1:N
				track = tracks(i);
	
				if (i == 1)
					tau = track.tau;
				elseif (tau ~= track.tau)
					disp('Error: tracks must have the same timestep');
					return;
				end

				X(:,i) = track.steps(:,1);
				Y(:,i) = track.steps(:,2);
				Z(:,i) = track.steps(:,3);
			end
		end

	end % END Static methods

	methods
		
		% constructor
		function [self] = HMMSolver(S, options) 
			self.S = S;
			if (nargin == 2)
				self.options = options;
			end
		end

		function [D, V] = solve(self, track)
			if (~ isempty(self.hmms))
				[D, V] = self.infer(track);
			else
				fprintf('Did you forget to train the solver?\n');
				return
			end
		end

		function [errs] = compare(self, track, D, V)
			errs = track.compare(D(:,1), V);
		end

		% Select the ML HMM that explains the observed tracks.
		% Require all tracks to have the same timestep.
		% 
		% @param tracks 	mat					1xN
		% @param maxIter	int					number of EM iterations
		%
		%	@return	LL				vec					log-likelihood in each dimension
		function LL = train(self, tracks)
			[X,Y,Z] = HMMSolver.readTracks(tracks);
			self.hmms = CHMM.empty(3,0);
			LL = zeros(3, 1);
			[self.hmms(1), LL(1)] = CHMM.fit(self.S, X, self.options);
			[self.hmms(2), LL(2)] = CHMM.fit(self.S, Y, self.options);
			[self.hmms(3), LL(3)] = CHMM.fit(self.S, Z, self.options);
		end

	end
end
