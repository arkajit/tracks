classdef TrackSolver < handle

	methods (Abstract)
		[D, V] = solve(self, track);
		[errs] = compare(self, track, D, V);
	end

	properties (Constant, Abstract)
		isState
	end

	methods
	
		function [res, errs] = test(self, tracks)
			N = length(tracks);
			errs.AE = nan(N, 4);
			errs.RE = nan(N, 4);	
			res(N).D = [];
			res(N).V = [];
			for i=1:N
				t = tracks(i);
				[res(i).D, res(i).V] = self.solve(t);
				errs.AE(i,:) = self.compare(t, res(i).D, res(i).V);
				if (self.isState)
					errs.RE(i,:) = t.compareStates(res(i).D(:,1), res(i).V);
				end
			end
		end 
	end

	methods (Static)
		function [results] = testAll(solvers, tracks)
			S = length(solvers);
			N = length(tracks);
			res(S, N).D = [];
			res(S, N).V = [];
			for i=1:S
				fprintf('Testing solver %d\n', i);
				s = solvers(i);
				res(i,:) = s.test(tracks);
			end
		end
	end	

end
