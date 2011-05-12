classdef TrackSolver < handle

	methods (Abstract)
		[D, V] = solve(self, track);
		[errs] = compare(self, track, D, V);
	end

	methods
		function [res] = test(self, tracks)
			N = length(tracks);
			for i=1:N
				t = tracks(i);
				[res(i).D, res(i).V] = self.solve(t);
				[res(i).errs] = self.compare(t, res(i).D, res(i).V);
			end
		end 
	end

	methods (Static)
		function solveAll(track, solvers)
			tp = TrackPlotter(track);
			for s=solvers
				[D, V] = s.solve(track);
				tp.plotDiffusion(D);
				tp.plotVelocity(V);
			end
		end
	end	

end
