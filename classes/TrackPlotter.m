classdef TrackPlotter

	properties
		track
	end
	
	methods

		function [self] = TrackPlotter(track)
			self.track = track;
		end

		function [] = plotDiffusion(self, D)
			X = 1:(self.track.T - 1);
			
			figure;
			plot(X, [self.track.D(X), D(X)]);
			legend('Actual', 'Estimated');
      ylabel('Diffusion Coefficient ([dist]^2/s)')
      xlabel('n (# steps)');
		end

		function [] = plotVelocity(self, V, useMag)
			if (nargin == 2)
				useMag = false;
			end;

			X = 1:(self.track.T - 1);

			if (useMag)
				figure;
				plot(X, [self.track.normV(X), V(X)]);
				legend('Actual', 'Estimated');
      	ylabel('Velocity Magnitude ([dist]/s)');
      	xlabel('n (# steps)');
				return;
			end

			labels = 'xyz';
			for i=1:3	
				figure;
				plot(X, [self.track.V(:,i), V(:,i)]);
				legend(sprintf('v_{%s} (Actual)', labels(i)), ...
							 sprintf('v^*_{%s} (Estimated)', labels(i)));
				ylabel(sprintf('Instantaneous v_{%s} ([dist]/s)', labels(i)));
				xlabel('n (# steps)');
			end
		end

	end

end
