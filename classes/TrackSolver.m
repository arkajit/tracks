classdef TrackSolver < handle

	methods (Abstract)
		[D, V] = solve(self, track);
		[errs] = compare(self, track, D, V);
	end

	methods
		function [res] = test(self, tracks)
			N = length(tracks);
			res(N).D = [];
			res(N).V = [];
			res(N).errs = [];
			for i=1:N
				t = tracks(i);
				[res(i).D, res(i).V] = self.solve(t);
				[res(i).errs] = self.compare(t, res(i).D, res(i).V);
			end
		end 
	end

	methods (Static)
		function [results] = testAll(solvers, tracks)
			S = length(solvers);
			N = length(tracks);
			res(S, N).D = [];
			res(S, N).V = [];
			res(S, N).errs = [];
			for i=1:S
				s = solvers(i);
				res(i,:) = s.test(tracks);
			end
		end
	end	

end
