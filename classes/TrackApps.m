classdef TrackApps
	properties (Constant)
		T = 100;
	end

	methods (Static)
    function app = diff()
			T = TrackApps.T;
			app.D = ones(T, 1);
			app.V = zeros(T, 3);
			app.tau = 1;
		end

    function app = diffstep(d)
			T = TrackApps.T;
			app.D = [ones(0.5*T, 1); d*ones(0.5*T, 1)];
			app.V = zeros(T, 3);
			app.tau = 1;
		end

		% normalize D to 1, xv = delta v_x
		function app = velstep(xv)
			T = TrackApps.T;
			app.D = ones(T, 1);
			app.V = [[zeros(0.5*T, 1); xv*ones(0.5*T, 1)] [zeros(T,2)]];
			app.tau = 1;
		end

    function app = diffvelstep(d, xv)
			T = TrackApps.T;
			app.D = [ones(0.25*T, 1); d*ones(0.75*T, 1)];
			app.V = [[zeros(0.6*T, 1); xv*ones(0.4*T, 1)] [zeros(T,2)]];
			app.tau = 1;
		end

    function app = pulse(d, xv, t)
			T = TrackApps.T;
			t1 = 0.3*T;
			t2 = t1+t;
			app.D = [ones(t1, 1); ones(t, 1); ones(T-t2, 1)];
			app.V = [[zeros(t1, 1); xv*ones(t,1); ones(T-t2, 1)] [zeros(T,2)]];
			app.tau = 1;
    end
	end
end
