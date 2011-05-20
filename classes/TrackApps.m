classdef TrackApps
	properties 
		T
		tau
	end

	methods
	
		function [self] = TrackApps(T, tau)
			if (nargin < 1)
				self.T = 100;
			else
				self.T = T;
			end

			if (nargin < 2)
				self.tau = 1;
			else
				self.tau = tau;
			end
		end
	
    function app = diff(self)
			T = self.T;
			app.name = 'DIFF';
			app.D = ones(T, 1);
			app.V = zeros(T, 3);
			app.tau = self.tau;
		end

    function app = diffstep(self, d)
			T = self.T;
			app.name = 'DIFFSTEP';
			app.D = [ones(0.5*T, 1); d*ones(0.5*T, 1)];
			app.V = zeros(T, 3);
			app.tau = self.tau;
		end

		% normalize D to 1, xv = delta v_x
		function app = velstep(self, xv)
			T = self.T;
			app.name = 'VELSTEP';
			app.D = ones(T, 1);
			app.V = [[zeros(0.5*T, 1); xv*ones(0.5*T, 1)] [zeros(T,2)]];
			app.tau = self.tau;
		end

    function app = diffvelstep(self, d, xv)
			T = self.T;
			app.name = 'DIFFVELSTEP';
			app.D = [ones(0.25*T, 1); d*ones(0.75*T, 1)];
			app.V = [[zeros(0.6*T, 1); xv*ones(0.4*T, 1)] [zeros(T,2)]];
			app.tau = self.tau;
		end

    function app = pulse(self, xv, t)
			T = self.T;
			t1 = 0.25*T;
			t2 = t1+t;
			app.name = 'PULSE';
			app.D = ones(T, 1);
			app.V = [[zeros(t1, 1); xv*ones(t,1); zeros(T-t2, 1)] [zeros(T,2)]];
			app.tau = self.tau;
    end
	end
end
