classdef TrackAnalyzer
  % TrackAnalyzer an object for running controlled experiments over several
  % tracks by varying different parameters

  properties (Constant=true)
    sdate='HH:MM_mmmdd';    % date format appended to output (11:26_Dec24)
  end

  properties
    k       % number of trials to average error measurements over
    nbins   % number of bins to use in track recovery
    ntrans  % number of times to transition in random test trajectories
  end

  methods
    % constructor
    function self = TrackAnalyzer(k, nbins, ntrans)
      self.k = k;
      self.nbins = nbins;
      self.ntrans = ntrans;
    end

    % analyze fixed V^2 *tau / 2D ratio as you vary track size
    % normalize D to 1
    function errS = analyzeLength(self, Smax, Vmax, tau)
      Dmax = 1;
      errS = zeros(Smax,1);

      for s=2*self.ntrans:Smax
        if mod(s,10) == 0
          disp(sprintf('s = %d', s)); % progress report
        end

        for i=1:self.k
          track = Track.fixedRandTrack(s, self.ntrans, Dmax, Vmax, tau);
          errS(s) = errS(s) + track.recover(Dmax, Vmax, self.nbins);
        end
      end

      errS = errS ./ self.k;
      save(sprintf('errS_%s.mat', datestr(now, TrackAnalyzer.sdate)), ...
                   'self', 'errS');
      plot(errS);
    end

    % analyze V, D, and S as you vary tau over taus
    % normalize D to 1
    function errTau = analyzeTau(self, S, Vmax, taus)
      Dmax = 1;
      errTau = zeros(length(taus), 1);

      for n=1:length(taus)
        if mod(n, 10) == 0
          disp(sprintf('n = %d', n)); % progress report
        end

        tau = taus(n);

        for i=1:self.k
          track = Track.fixedRandTrack(S, self.ntrans, Dmax, Vmax, tau);
          errTau(n) = errTau(n) + track.recover(Dmax, Vmax, self.nbins);
        end
      end

      errTau = errTau ./ self.k;
      save(sprintf('errTau_%s.mat', datestr(now, TrackAnalyzer.sdate)), ...
           'self', 'errTau', 'taus');
      plot(taus, errTau);
    end

    % analyze fixed S, D, and tau as you vary Vmax over vels
    % normalize D to 1
    function errRat = analyzeRatio(self, S, vels, tau)
      Dmax = 1;
      errRat = zeros(length(vels),1);
     
      for n=1:length(vels)
        if mod(n, 10) == 0
          disp(sprintf('n = %d', n)); % progress report
        end

        Vmax = vels(n);

        for i=1:self.k
          track = Track.fixedRandTrack(S, self.ntrans, Dmax, Vmax, tau);
          errRat(n) = errRat(n) + track.recover(Dmax, Vmax, self.nbins);
        end
      end

      errRat = errRat ./ self.k;
      save(sprintf('errRat_%s.mat', datestr(now, TrackAnalyzer.sdate)), ...
           'self', 'errRat', 'vels');
      plot(vels, errRat);
    end

    function deltaTs = analyzeOneTrans(self, deltaVs)
      D = 1;
      tau = 1;
      S = 50;
      deltaTs = zeros(length(deltaVs), 1);
    
      for n=1:length(deltaVs)
       
        deltaV = deltaVs(n); 
        vel = [deltaV 0 0]; % just turn on velocity in x direction

        for i=1:self.k
          t = Track.randomTrack(D*ones(S, 1), zeros(S,3), tau);
          s = Track.randomTrack(D*ones(S, 1), repmat(vel, S, 1), tau);
          r = t+s;
          [~, Ds, Vs] = r.recover(D, deltaV, self.nbins);

          vx = Vs(:,1); % just care about the x
          [b, m] = unique(vx);

          if (size(b, 1) == 2)
            deltaT = tau * abs(m(1) - S); % S is the true last index of old regime
          elseif (size(b,1) == 1)
            deltaT = tau*S; % transition not found
            % penalize with highest absolute error
          else
            deltaT = tau*S; % too many transitions found
            disp(sprintf('Not Two: %d', size(b,1)));
            disp(b);
            disp(m);
            % TODO(arkajit): fix this -- trying to find the segment that's
            % closest to the true value and then discover where the previous
            % segment ended OR just find a better way
            %[~, idx] = min(abs(b - deltaV));
            %m = sort(m);
            %if idx > 1
            %  deltaT = tau * abs(m(idx-1) - S);
          end

          deltaTs(n) = deltaTs(n) + deltaT;

        end 
      end

      deltaTs = deltaTs ./ self.k;
    end

  end
end
