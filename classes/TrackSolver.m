classdef TrackSolver
	
	methods (Abstract)
		[D, V] = solve(self, track);	
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
