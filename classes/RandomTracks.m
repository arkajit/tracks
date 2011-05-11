classdef RandomTracks
  % RandomTracks is a utility class for generating random tracks.

  methods (Static, Access = private)

    % Generate a random trajectory from sequence of motion parameters D and V. 
    % The length of the sequences should be equal and will be the number of
    % steps generated in the trajectory. If the number of steps in the
    % trajectory is denoted S:
    % 
    % @param D	 	mat 		Sx1
    % @param V 		mat 		Sx3
		% @param tau	double 	the time between steps
    %
    % @return pos mat (S+1)x3 
    function [pos] = generateTrajectory(D, V, tau)
      S = size(D,1); % number of steps
      T = S+1;       % number of points

      if size(V,1) ~= S 
        pos = []; % dimensions don't match, abort
        return;
      else
        pos = zeros(T,3);
        pos(1,:) = rand(1,3); % start at random point in unit box
      end

      for t=1:S
        pos(t+1,:) = pos(t,:) + sqrt(2*D(t)*tau)*randn(1,3) + tau*V(t,:);
      end
    end

  end

  methods (Static)
   
    % Generates a random @Track from a sequence of motion parameters, D, V.
    %
    % @param D mat
    % @param V mat
    % @param tau double  
    %
    % @return track @Track
    function [track] = from(D, V, tau)
      trajectory = RandomTracks.generateTrajectory(D, V, tau);
      track = Track.fromTrajectory(trajectory, tau);
      track.D = D;
      track.V = V;
    end

		function [track] = random(nSteps)
			if (nargin == 0)
				nSteps = 100;
			end
			Dmax = randi(20) / 10;
			Vmax = randi(100) / 10;
			tau = 1;
			track = RandomTracks.fromParams(nSteps, Dmax, Vmax, tau);
		end

		function [tracks] = sample(hmms, tau, N)
			if (nargin < 3)
				N = 100;
			end

			for i=1:N
				tracks(i) = RandomTracks.fromModels(hmms, tau);
			end
		end

		function [track] = fromModels(hmms, tau, T)
			if (nargin < 3)
				T = 100;
			end

			[~, Sx] = hmms(1).sample(1, T);
			[~, Sy] = hmms(2).sample(1, T);
			[~, Sz] = hmms(3).sample(1, T);

			D = zeros(T, 1);
			D = (hmms(1).stddevs(Sx) .^ 2) ./ (2 * tau);
			
			V = zeros(T, 3);
			V(:,1) = hmms(1).means(Sx) ./ tau;
			V(:,2) = hmms(2).means(Sy) ./ tau;
			V(:,3) = hmms(3).means(Sz) ./ tau;
			track = RandomTracks.from(D, V, tau);
		end

    function [track] = diffusion(nSteps, diffusionCoeff, tau)
      D = diffusionCoeff*ones(nSteps,1);
      V = zeros(nSteps, 3);
      track = RandomTracks.from(D, V, tau);
    end
    
    function [track] = flow(nSteps, velocity, tau)
      D = zeros(nSteps, 1);
      V = repmat(velocity, nSteps, 1);
      track = RandomTracks.from(D, V, tau); 
    end

    function [track] = diffStep(nSteps, d1, d2, tau) 
      n = randi([floor(nSteps/4) floor(3*nSteps/4)]); % random step point

      D1 = d1*ones(n,1);
      V1 = zeros(n, 3);
      t1 = RandomTracks.from(D1, V1, tau);

      D2 = d2*ones(nSteps-n,1);
      V2 = zeros(nSteps-n, 3);
      t2 = RandomTracks.from(D2, V2, tau);
    
      track = t1+t2;
    end

    function [track] = step(nSteps, diffusionCoeff, xvel1, xvel2, tau) 
      n = randi([floor(nSteps/4) floor(3*nSteps/4)]); % random step point

      D1 = diffusionCoeff*ones(n,1);
      V1 = repmat([xvel1 0.1 0.1], n, 1);
      t1 = RandomTracks.from(D1, V1, tau);

      D2 = diffusionCoeff*ones(nSteps-n,1);
      V2 = repmat([xvel2 0.1 0.1], nSteps-n, 1);
      t2 = RandomTracks.from(D2, V2, tau);
    
      track = t1+t2;
    end

    function [track] = fromParams(nSteps, Dmax, Vmax, tau, smooth) 
      D = Dmax*rand(nSteps,1);
      V = 2*Vmax*rand(nSteps,3)-Vmax;

			if (nargin < 5)
				smooth = true;
			else
				smooth = false;
			end
		
			if (smooth)	
				% smooth out transitions
				i = 1;
				while i < 0.9*nSteps
					j = randi([i,nSteps]);
					D(i:j-1) = D(i);
					i = j;
				end
				D(i:end) = D(i);

				% do V transitions independently of D
				i = 1;
				while i < 0.9*nSteps
					j = randi([i,nSteps]);
					V(i:j-1,:) = repmat(V(i,:), j-i, 1);
					i = j;
				end
				V(i:end,:) = repmat(V(i,:), nSteps-i+1, 1);
			end

      track = RandomTracks.from(D, V, tau);
    end 

    function [track] = fixedTransitions(nSteps, nTrans, Dmax, Vmax, tau)
      D = Dmax*rand(nSteps,1);
      V = 2*Vmax*rand(nSteps,3)-Vmax;

      r = randperm(nSteps);
      nums = sort(r(1:nTrans));
      last = 1;
      for i=nums
        D(last:i-1) = D(last);
        last = i;
      end;
      D(last:end) = D(i);
   
      r = randperm(nSteps);
      nums = sort(r(1:nTrans));
      last = 1;
      for i=nums
        V(last:i-1,:) = repmat(V(last,:), i-last, 1);
        last = i;
      end;
      V(last:end,:) = repmat(V(last,:), nSteps-last+1, 1);

      track = RandomTracks.from(D, V, tau);
    end

  end
end
