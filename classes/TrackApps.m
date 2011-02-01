classdef TrackApps

  properties (Constant)
    opts = statset('FunValCheck', 'off'); 
  end

  methods (Static)

    function [Ds, Vs] = testOneTrans(S, D, deltaV, tau)
      r = TrackApps.getOneTransTrack(S, D, deltaV, tau);
      [~, Ds, Vs]  = r.recover(D, deltaV, 20);
    end

    function [r] = getOneTransTrack(S, D, deltaV, tau)
      t = Track.randomTrack(D*ones(S,1), zeros(S,3), tau);
      s = Track.randomTrack(D*ones(S,1), repmat([deltaV 0 0], S, 1), tau);
      r = t+s;
    end

    % Implements Arcizet et. al.'s TRAnSpORT algorithm which computes a
    % rolling-average MSD over M frames.
    % 
    % @param track Track
    % @param M int - size of window in frames
    % @param sig_alpha double - tolerance on alpha to delcare active state
    % @param sig_phi double - tolerance on angle correlation to declare active
    % 
    % @return D - predicted diffusion coefficients
    % @return V - predicted velocity magnitudes
    % @return pA - predicted active/passive states
    % 
    % If the last three parameters are not supplied, will use the defaults used
    % in Arcizet et. al.
    % Instead of computing MSD versus time (t), do it versus step # (n).
    % NOTE(arkajit): Thus this implicitly assumes that tau = 1       
    function [D, V, pA] = TRANSPORT(track, M, sig_alpha, sig_phi)
      if (nargin == 1)
        M = 40;
        sig_alpha = 0.3;
        sig_phi = 0.6;
      end

      R = track.positions;
      T = size(R, 1);

      D = zeros(T, 1);
      V = zeros(T, 1);
      pA = zeros(T, 1);
    
      for i=1:T
        beta = TrackApps.plfit(R, i, M, [1; 1]);  % arbitrary initial condition
        A = beta(1);
        alpha = beta(2);
        phi = TrackApps.anglecorr(track.angles, i, M);
        if (all([alpha >= 2 - sig_alpha, alpha <= 2 + sig_alpha, ...
                 phi >= -sig_phi, phi <= sig_phi]))
          pA(i) = 1;
          V(i) = sqrt(A);
        else
          D(i) = A/6; % 6 for 3D and assumes tau = 1
        end
      end
    end

    % For a trajectory matrix, do a polynomial fit of its MSD curve over a
    % certain subwindow.
    function [beta] = plfit(R, i, M, b0)
      msdR = TrackApps.msds3D(R, i, M);
      X = (1:size(msdR, 1))';
      beta = nlinfit(X, msdR, @powerlaw, b0, TrackApps.opts);
      plot(X, [msdR, powerlaw(beta, X)]);

      function [yhat] = powerlaw(b, X)
        A = b(1);
        alpha = b(2);
        yhat = A * X.^alpha;
      end
    end

    function [stddev] = anglecorr(angles, i, M)
      T = size(angles, 1);
      lower = max(1, i - floor(M/2));
      upper = min(T, i + floor(M/2));

      msdFunc = TrackApps.getMSD(angles(lower:upper));
      stddev = msdFunc(floor(M/4)); % Arcizet et. al. use M/4 as their lag
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
      msds = TrackApps.allMSDs3D(R(lower:upper,:));

      % 3D MSD is sum of MSDs in each direction since random walkers are independent
      msdR = sum(msds, 2); 
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
        msds(:,i) = TrackApps.allMSDs(R(:,i));
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
      msdFunc = TrackApps.getMSD(x);
      for i=1:T-1
        msds(i) = msdFunc(i);
      end
    end

    % For any set of data points x, compute and return the MSD function for it.
    % This computes a time average within an individual trajectory, NOT an
    % ensemble average over a set of trajectories.
    % 
    % @param x - trajectory vector
    % @return msd - function mapping number of steps, n, to the mean square
    %               displacement of x over n steps
    function [msd] = getMSD(x)
      T = length(x);
      msd = @my_msd;     
 
      function [s] = my_msd(n)
        vals = zeros(T-n, 1);
        for i=1:T-n
          vals(i) = (x(i+n) - x(i))^2;
        end
        s = mean(vals);
      end
    end

    function [nD, nV] = numStates(D, V)
      nD = length(unique(D));
      nV = [length(unique(V(:,1))); length(unique(V(:,2))); ...
            length(unique(V(:,3)))];
    end

    function [deltaN] = getDelta(V1, V2)
      [~, nV1] = TrackApps.numStates(zeros(1), V1);
      [~, nV2] = TrackApps.numStates(zeros(1), V2);
      if ((nV1(1) == nV2(1)) && (nV1(1) == 2))
        [~, m1] = unique(V1(:,1));
        [~, m2] = unique(V2(:,1));
        deltaN = abs(m2(1) - m1(1)); 
      else
        deltaN = NaN;
      end
    end

  end

end
