classdef HMMSolver < handle
	% represents the logic used to solve a Track using an HMM

	properties
		%Dmax
		%Vmax
		S 		% expected number of states
		hmm
		hmms
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
					hmm = self.hmms{i};
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
      %maxstd = sqrt(2*self.Dmax*tau);
      %maxmean = self.Vmax*tau;
			hmm = CHMM.random(self.S);
		end

	end

	methods (Static)

		function [frac] = pctErr(A, B)
			frac = 100 * norm(A-B) / norm(A);
		end

		function [frac, errs, order] = errors(n, A, B)
			errs = -1;
			alen = length(A);
			order = [];
			for f=perms(1:n)'
				nErrs = 0;
				
				% compute the errors with this permutation
				for i=1:alen
					nErrs = nErrs + sum(A{i} ~= f(B{i}));
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

			frac = errs / (alen * length(A{1}));
		end

	end

	methods
		
		% constructor
		function [self] = HMMSolver(S) 
			%self.Dmax = Dmax;
			%self.Vmax = Vmax;
			self.S = S;
		end

		% Select the ML HMM that explains the observed tracks.
		% Require all tracks to have the same timestep.
		% 
		% @param tracks 	cellarray 	Nx1
		% @param maxIter	int					number of EM iterations
		function train(self, tracks, maxIter)
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

		function trainXYZ(self, tracks)
			N = length(tracks);
			if (~N)
				return;
			end

			X = cell(N, 1);
			Y = cell(N, 1);
			Z = cell(N, 1);

			for i=1:N
				track = tracks{i};
	
				if (i == 1)
					tau = track.tau;
				elseif (tau ~= track.tau)
					disp('Error: tracks must have the same timestep');
					return;
				end

				X{i} = track.steps(:,1);
				Y{i} = track.steps(:,2);
				Z{i} = track.steps(:,3);
			end

			self.hmms = cell(3,1);
			self.hmms{1} = CHMM.fit(self.S, X);
			self.hmms{2} = CHMM.fit(self.S, Y);
			self.hmms{3} = CHMM.fit(self.S, Z);
		end

		function [D, V, errs] = test(self, tracks)
			N = length(tracks);		% no. of test examples
			D = cell(N, 1);
			V = cell(N, 1);
			errs = zeros(N, 4);

			for i=1:N
				t = tracks{i};
				[D{i}, V{i}] = self.infer(t);
				errs(i,1) = HMMSolver.pctErr(t.D, D{i}(:,1));
				for j=1:3
					errs(i,j+1) = HMMSolver.pctErr(t.V(:,j), V{i}(:,j));
				end
			end
		end

	end
end
