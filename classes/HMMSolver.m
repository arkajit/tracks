classdef HMMSolver < TrackSolver
	% represents the logic used to solve a Track using an HMM

	properties
		Dmax
		Vmax
		P			% expected number of velocity states
		Q			% expected number of diffusion states
		hmm
	end	

	methods (Access = private)
	
		% Using selected HMM, perform inference to find the most likely sequence of
		% motion parameters that could generate this track. Computes a Viterbi
		% decoding.
		% 
		% @param track Track
		% 
		% @return D vec
		% @return V mat 
		function [D3, V] = infer(self, track)
      D = zeros(size(track.D));
      D3 = zeros(size(track.V)); % want 3D
      V = zeros(size(track.V));

      for i=1:3
				states = self.hmm.viterbi(track.steps(:,i));

				V(:,i) = hmm.means(states) ./ track.tau;
				D3(:,i) = (hmm.stddevs(states) .^ 2) ./ (2 * track.tau);
      end
		
      %D = mean(D3, 2); 	% average x,y,z predictions
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
		function [hmm] = initFHMM(self, tau)
      maxstd = sqrt(2*self.Dmax*tau);
      maxmean = self.Vmax*tau;
			odds = 10;					% how to initialize the transition matrix

			sigs = rand(self.Q, 1) * maxstd;	
			mus = rand(self.P, 1) .* sign(randn(self.P, 1)) * maxmean;

			S = self.P*self.Q;	% number of states
			pi = ones(S, 1);
			if S == 1
				M = [1];
			else
				M = ones(S) + ((S-1)*odds - 1)*eye(S); % diag. entries = odds*(S-1)
			end
			
			mc = MarkovChain(pi, M);
			means = reshape(repmat(mus, Q, 1), S, 1);
			stddevs = reshape(repmat(sigs', P, 1), S, 1);
			hmm = CHMM(mc, means, stddevs);
		end

	end

	methods
		
		% constructor
		function [self] = HMMSolver(Dmax, Vmax, Q, P) 
			self.Dmax = Dmax;
			self.Vmax = Vmax;
			self.Q = Q;
			self.P = P;
		end

		% Select the ML HMM that explains the observed tracks.
		% Require all tracks to have the same timestep.
		% 
		% @param tracks 	cellarray 	Nx1
		function train(self, tracks)
			N = length(tracks);
			if (!N)
				return;
			end

			X = cell(3*N, 1);
			tau = 0;
			for i=1:N
				track = tracks{i};

				if (!tau) 
					tau = track.tau;
				else if (tau != track.tau)
					disp('Error: tracks must have the same timestep');
					return;
				end

				for j=1:3
					X{i*j} = track.steps(:,j);
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

	end
end
