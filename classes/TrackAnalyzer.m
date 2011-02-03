classdef TrackAnalyzer
  % TrackAnalyzer an object for running controlled experiments over several
  % tracks by varying different parameters

  properties (Constant=true)
    sdate='HH:MM_mmmdd';    % date format appended to output (11:26_Dec24)
  end

  properties
    nTrials   % number of trials to average error measurements over
    nBins     % number of bins to use in track recovery
    nTrans    % number of times to transition in random test trajectories
  end

  methods

    function self = TrackAnalyzer(nTrials, nBins, nTrans)
      self.nTrials = nTrials;
      self.nBins = nBins;
      self.nTrans = nTrans;
    end

    % Analyze fixed V^2 * tau / 2D ratio while varying track length. Normalizes
    % diffusion coefficient to 1. Returns data to plot error versus track
    % length.
    %
    % @param Smax int - max # of steps in track
    % @param Vmax double - max velocity
    %
    % @return errS vec (Smax)x1
    function errS = analyzeLength(self, Smax, Vmax, tau)
      Dmax = 1;
      errS = zeros(Smax,1);

      for s=2*self.nTrans:Smax
        if mod(s,10) == 0
          disp(sprintf('s = %d', s)); % progress report
        end

        for i=1:self.nTrials
          track = RandomTracks.fixedTransitions(s, self.nTrans, Dmax, Vmax, tau);
          errS(s) = errS(s) + track.recover(Dmax, Vmax, self.nBins);
        end
      end

      errS = errS ./ self.nTrials;
      save(sprintf('errS_%s.mat', datestr(now, TrackAnalyzer.sdate)), ...
                   'self', 'errS');
      plot(errS);
    end

    % Analyze V, D, and S as you vary tau over a set of taus. Diffusion
    % coefficient is normalized to 1. Returns data to plot error versus time
    % lag, tau (holding other factors constant).
    % 
    % @param S int - length of random track to use for analysis
    % @param Vmax double - max velocity
    % @param taus vec Lx1 - list of taus to try
    % 
    % @return errTau vec Lx1
    function errTau = analyzeTau(self, S, Vmax, taus)
      Dmax = 1;
      errTau = zeros(length(taus), 1);

      for n=1:length(taus)
        if mod(n, 10) == 0
          disp(sprintf('n = %d', n)); % progress report
        end

        tau = taus(n);

        for i=1:self.nTrials
          track = RandomTracks.fixedTransitions(S, self.nTrans, Dmax, Vmax, tau);
          errTau(n) = errTau(n) + track.recover(Dmax, Vmax, self.nBins);
        end
      end

      errTau = errTau ./ self.nTrials;
      save(sprintf('errTau_%s.mat', datestr(now, TrackAnalyzer.sdate)), ...
           'self', 'errTau', 'taus');
      plot(taus, errTau);
    end

    % Analyze fixed S, D, and tau as you vary Vmax over a set of velocities.
    % Diffusion coefficient is normalized to 1. Returns data to plot error
    % versus the maximum velocity Vmax.
    % 
    % @param S int - length of random track to use for analysis
    % @param vels vec Lx1 - list of max velocities to try
    %
    % @return errRat vec Lx1
    function errRat = analyzeRatio(self, S, vels, tau)
      Dmax = 1;
      errRat = zeros(length(vels),1);
     
      for n=1:length(vels)
        if mod(n, 10) == 0
          disp(sprintf('n = %d', n)); % progress report
        end

        Vmax = vels(n);

        for i=1:self.nTrials
          track = RandomTracks.fixedTransitions(S, self.nTrans, Dmax, Vmax, tau);
          errRat(n) = errRat(n) + track.recover(Dmax, Vmax, self.nBins);
        end
      end

      errRat = errRat ./ self.nTrials;
      save(sprintf('errRat_%s.mat', datestr(now, TrackAnalyzer.sdate)), ...
           'self', 'errRat', 'vels');
      plot(vels, errRat);
    end

  end % end public methods
end
