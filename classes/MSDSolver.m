classdef MSDSolver < TrackSolver

	properties (Constant)
    opts = statset('FunValCheck', 'off'); 
	end

	properties
		isRolling
		M						% window size in frames
		sigAlpha
		sigPhi	
	end

	methods (Static)

    % For any set of data points x, compute and return the MSD function for it.
    % This computes a time average within an individual trajectory, NOT an
    % ensemble average over a set of trajectories.
    % 
    % @param x - trajectory vector
    % @return msd - function mapping number of steps, n, to the mean square
    %               displacement of x over n steps
    function [msd] = getMSD(x)
      T = length(x);
      msd = @myMSD;     
 
      function [s] = myMSD(n)
        vals = zeros(T-n, 1);
        for i=1:T-n
          vals(i) = (x(i+n) - x(i))^2;
        end
        s = mean(vals);
      end
    end

    % For a given trajectory of T points, compute all MSDs for points that are
    % 1, 2, ..., T-1 steps apart within the trajectory.
    %
    % @param x - trajectory vector, Tx1
    % @return msds - vector of msds, (T-1)x1
    function [msds] = allMSDs(x)
      T = length(x);
      msds = zeros(T-1,1);
      msdFunc = MSDSolver.getMSD(x);
      for i=1:T-1
        msds(i) = msdFunc(i);
      end
    end

    % For a 3D trajectory, consider the x,y, and z trajectories independently
    % and compute all the MSDs for each.
    %
    % @param R mat - 3D trajectory matrix, Tx3 with each column representing a
    %            separate trajectory
    % @return msds - matrix of msds, (T-1)x3
    function [msds] = allMSDs3D(R)
      T = length(R);
      msds = zeros(T-1,3);
      for i=1:3
        msds(:,i) = MSDSolver.allMSDs(R(:,i));
      end
    end

    % For a 3D trajectory, compute the MSDs for a particular window of points
    % within it.
    %
    % @param R mat - 3D trajectory matrix
    % @param i int - frame index for the window's center
    % @param M int - length of the window in frames
    % @return msdR vec - 3D MSDs for each possible step size
    function [msdR] = msds3D(R, i, M)
      T = size(R, 1);
      lower = max(1, i - floor(M/2));
      upper = min(T, i + floor(M/2));
      msds = MSDSolver.allMSDs3D(R(lower:upper,:));

      % 3D MSD is sum of MSDs in each direction since random walkers are independent
      msdR = sum(msds, 2); 
    end

    % For a trajectory matrix, do a polynomial fit of its MSD curve over a
    % certain subwindow.
    function [beta] = plfit(R, i, M, fitFunc)
			% start with an arbitrary initial condition
			b0 = [1; 1];

      msdR = MSDSolver.msds3D(R, i, M);
      X = (1:size(msdR, 1))';
      beta = nlinfit(X, msdR, fitFunc, b0, MSDSolver.opts);
      plot(X, [msdR, fitFunc(beta, X)]);
    end

    function [stddev] = angleCorrelation(angles, i, M)
      T = size(angles, 1);
      lower = max(1, i - floor(M/2));
      upper = min(T, i + floor(M/2));

      msdFunc = MSDSolver.getMSD(angles(lower:upper));

			% Arcizet et. al. use M/4 as their lag
      stddev = msdFunc(floor(M/4)); 
    end

	end

	methods (Access = private)

    % Implements Arcizet et. al.'s TRAnSpORT algorithm which computes a
    % rolling-average MSD over M frames.
    % 
    % @param track Track
    % 
    % @return D - predicted diffusion coefficients
    % @return V - predicted velocity magnitudes
    % @return pA - predicted active/passive states
    % 
    % Instead of computing MSD versus time (t), do it versus step # (n).
    % NOTE(arkajit): Thus this implicitly assumes that tau = 1       
		function [D, V, pA] = TRANSPORT(self, track)
      R = track.positions;
      T = size(R, 1);

      D = zeros(T, 1);
      V = zeros(T, 1);
      pA = zeros(T, 1);

			function [yhat] = fitFunc(b, X)
				A = b(1);
				alpha = b(2);
				yhat = A * X.^alpha;
			end
    
      for i=1:T
        beta = MSDSolver.plfit(R, i, self.M, @fitFunc);
        A = beta(1);
        alpha = beta(2);

        phi = MSDSolver.angleCorrelation(track.angles, i, self.M);
        if (all([alpha >= 2 - self.sigAlpha, ...
								 alpha <= 2 + self.sigAlpha, ...
                 phi >= -self.sigPhi, ...
								 phi <= self.sigPhi
								]))
          pA(i) = 1;
          V(i) = sqrt(A);
        else
          D(i) = A/6; % 6 for 3D and assumes tau = 1
        end
      end
		end

		function [D, V] = MSD(self, track)
      R = track.positions;
      T = size(R, 1);

      D = zeros(T, 1);
      V = zeros(T, 1);

			function [yhat] = fitFunc(b, X)
				A = b(1);
				B = b(2);
				yhat = A * X + B * X.^2;
			end
  
      for i=1:T
        beta = MSDSolver.plfit(R, i, self.M, @fitFunc);
        A = beta(1);
        B = beta(2);
        D(i) = A/6;
        V(i) = sqrt(B);
      end

		end

	end

	methods

		% @param isRolling bool - whether to use the rolling MSD (TRANSPORT)
		% 												algorithm or the regular MSD
    % @param M int - size of window in frames
    % @param sigAlpha double - tolerance on alpha to delcare active state
    % @param sigPhi double - tolerance on angle correlation to declare active
		% If the last 3 parameters are not provided, we use defaults from Arcizet
		% et. al.
		function [self] = MSDSolver(isRolling, M, sigAlpha, sigPhi)
			if (nargin <= 2) 
				self.sigAlpha = 0.3;
				self.sigPhi = 0.6;
			else 
				self.sigAlpha = sigAlpha;	
				self.sigPhi = sigPhi;
			end	

			if (nargin <= 1)
				self.M = 40;
			else
				self.M = M;
			end

			if (nargin == 0)
				self.isRolling = false;
			else
				self.isRolling = isRolling;
			end
		end

		function [D, V] = solve(self, track)
			if (self.isRolling)
				[D, V] = self.TRANSPORT(track);
			else
				[D, V] = self.MSD(track);
			end		
		end

	end
end
