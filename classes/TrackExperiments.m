classdef TrackExperiments

  properties (Constant)
    S = 1000;
    nBins = 15;
    tau = 1;
  end

  methods (Static)

    function [t, D1, D2, V1, V2] = pureFlow()
      r = RandomTracks.flow(TrackExperiments.S/2, [-0.58 1.4 -1.76], ...
                            TrackExperiments.tau);
      s = RandomTracks.flow(TrackExperiments.S/2, [1.5, -0.7, 0.2], ...
                            TrackExperiments.tau);
      t = r+s;

      [D1, V1, D2, V2] = TrackExperiments.runAlgs(t, 0, 2);
      TrackExperiments.plotFlow(t, V1, V2);
    end

    function [t, D1, D2, V1, V2] = pureDiffusion()
      t = RandomTracks.diffusion(TrackExperiments.S, 1, TrackExperiments.tau);

      [D1, V1, D2, V2] = TrackExperiments.runAlgs(t, 1, 0);
      TrackExperiments.plotDiffusion(t, D1, D2);
    end

    function [t, D1, D2, V1, V2] = stepTrack()
      t = RandomTracks.step(TrackExperiments.S, 0.5, 0, 3, ...
                            TrackExperiments.tau);

      [D1, V1, D2, V2] = TrackExperiments.runAlgs(t, 1, 5);
      TrackExperiments.plotDiffusion(t, D1, D2);
      TrackExperiments.plotFlow(t, V1, V2);
    end

    function [t, D1, D2, V1, V2] = complexTrack(nTrans)
      t = RandomTracks.fixedTransitions(TrackExperiments.S, nTrans, 2, 5, ...
                                        TrackExperiments.tau);

      [D1, V1, D2, V2] = TrackExperiments.runAlgs(t, 2, 5);
      TrackExperiments.plotDiffusion(t, D1, D2);
      TrackExperiments.plotFlow(t, V1, V2);
    end

    function [D1, V1, D2, V2] = runAlgs(track, Dmax, Vmax)
      disp('Running MSD');
			m = MSDSolver();
      tic; [D2, V2] = m.solve(track); toc
      disp('Running HMM');
			h = HMMSolver(Dmax, Vmax, TrackExperiments.nBins);
      tic; [D1, V1] = h.solve(track); toc
    end

    function [] = plotDiffusion(track, D1, D2)
      X = 1:TrackExperiments.S;

      figure;
      plot(X, [track.D(X), D1(X), D2(X)]);
      legend('Actual', 'HMM', 'Rolling MSD');
      ylabel('Diffusion Coefficient ([dist]^2/s)')
      xlabel('n (# steps)');
    end

    function [] = plotFlow(track, V1, V2)
      X = 1:TrackExperiments.S;

      figure;
      plot(X, [V1, track.V]); 
      legend('v^*_x (HMM)', 'v^*_y (HMM)', 'v^*_z (HMM)',...
             'v_x (Actual)', 'v_y (Actual)', 'v_z (Actual)');
      ylabel('Instantaneous Velocities in each direction ([dist]/s)');
      xlabel('n (# steps)');

      figure;
      plot(X, [V2(X), track.normV(X)]);
      legend('Rolling MSD', 'Actual');
      ylabel('Velocity Magnitude ([dist]/s)');
      xlabel('n (# steps)');
    end

  end
end
