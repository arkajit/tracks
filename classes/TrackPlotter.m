classdef TrackPlotter

	properties
		track
		savedir 	% if set, will save plots here instead of showing them
	end
	
	methods

		function [self] = TrackPlotter(track, savedir)
			self.track = track;
			if nargin == 2
				if savedir(length(savedir)) ~= '/'
					self.savedir = [savedir '/'];
				else
					self.savedir = savedir;
				end
			else
				self.savedir = 0; % don't save, just show figure
			end
		end

		function [] = plotDiffusion(self, D)
			X = 1:(self.track.T - 1);
		
			if self.savedir	
				h = figure('visible', 'off');
			else
				h = figure;
			end
			plot(X, [self.track.D(X), D(X)]);
			legend('Actual', 'Estimated');
      ylabel('Diffusion Coefficient ([dist]^2/s)')
      xlabel('n (# steps)');

			if self.savedir
				saveas(h, [self.savedir 'diff.png']);
			end
		end

		function [] = plotVelocity(self, V, useMag)
			if (nargin == 2)
				useMag = false;
			end;

			X = 1:(self.track.T - 1);

			if (useMag)
				if (self.savedir)
					h = figure('visible', 'off');
				else
					h = figure;
				end
				plot(X, [self.track.normV(X), V(X)]);
				legend('Actual', 'Estimated');
      	ylabel('Velocity Magnitude ([dist]/s)');
      	xlabel('n (# steps)');
				
				if (self.savedir)
					saveas(h, [self.savedir 'vmag.png']);
				end
				return;
			end

			labels = 'xyz';
			for i=1:3
				if (self.savedir)	
					h = figure('visible', 'off');
				else
					h = figure;
				end;
				plot(X, [self.track.V(:,i), V(:,i)]);
				legend(sprintf('v_{%s} (Actual)', labels(i)), ...
							 sprintf('v^*_{%s} (Estimated)', labels(i)));
				ylabel(sprintf('Instantaneous v_{%s} ([dist]/s)', labels(i)));
				xlabel('n (# steps)');
			
				if (self.savedir)
					saveas(h, [self.savedir sprintf('v%s.png', labels(i))]);
				end
			end
		end

	end

end
