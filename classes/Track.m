classdef Track
  % Track represents an observed trajectory of a particle.
  %   Contains the steps and angles that describe the trajectory. Also may
  %   contain the true motion parameters, D and V, that generated the
  %   trajectory.
    
  properties
    T          % number of position observations
    positions  % T x 3 matrix of x, y, z positions
    steps      % (T-1) x 3 matrix of steps in x,y,z direction
    angles     % (T-2) x 1 matrix of angles between each step

    tau        % size of discrete timesteps (seconds)
    D          % (T-1) * 1 matrix of diffusion coeffs. (if known) for each step
    V          % (T-1) * 3 matrix of x,y,z velocities (if known) for each step
               % if D or V not known, will be a matrix of NaNs
  end
  
  properties (Dependent)
    normV
  end

  properties (Constant)
    epsilon = 1e-12;  % precision
    incr    = 1e4;    % resolution
  end

  methods (Static)

    function [theta] = angle(u,v)
      theta = acos(dot(u,v) / (norm(u) * norm(v)));
    end;

    % Round x only if it is within eps=1e-12 of the next floating point
    % where next is in increments of 1/incr
    function [y] = epsround(x)
      i = round(Track.incr*x)/Track.incr; % round to log(incr) decimal points
      if (abs(x - i) <= Track.epsilon)
        y = i;
      else
        y = x;
      end
    end

		function [frac] = percentError(A, B)
			frac = 100 * norm(A-B) / norm(A);
		end

    % From a trajectory (position time series), compute its length and the
    % sequence of steps and angles.
    % 
    % @param trajectory mat Tx3
    % 
    % @return T int - number of points in the trajectory
    % @return steps mat (T-1)x3 - steps in the trajectory
    % @return angles mat (T-2)x3 - angles in the trajectory
    function [T, steps, angles] = observe(trajectory)
      T = size(trajectory,1);
      steps = zeros(T-1, 3);
      angles = zeros(T-2, 1);

      for i=1:T-1
        % to prevent precision errors, round within an epsilon
        steps(i,:) = Track.epsround(trajectory(i+1,:) - trajectory(i,:));
      end

      for i=1:T-2
        u = steps(i,:);
        v = steps(i+1,:);
        angles(i) = Track.angle(u,v);
      end
    end

    function [self] = fromTrajectory(trajectory, tau)
      self = Track();
      self.positions = trajectory;
      self.tau = tau;
      [self.T, self.steps, self.angles] = Track.observe(trajectory);
      self.D = NaN(self.T-1,1);
      self.V = NaN(self.T-1,3);
    end

		% NOTE(arkajit): careful in calling this, steps should be drawn from a D, V
		% distribution, not just a mean, sigma distribution
		function [self] = fromSteps(steps, tau)
			self = Track();
			self.steps = steps;
			S = length(steps);
			self.T = S+1;
			pos = zeros(self.T, 3);

			for t=1:S
				pos(t+1,:) = pos(t,:) + steps(t,:);
			end

			self.positions = pos;
			[~, ~, self.angles] = Track.observe(pos);
			self.tau = tau;
      self.D = NaN(self.T-1,1);
      self.V = NaN(self.T-1,3);
		end

  end
    
  methods

    % empty constructor
    function self = Track()
      self.positions = [];
    end

    function nV = get.normV(self)
      nV = zeros(self.T-1, 1);
      for i=1:self.T-1
        nV(i) = norm(self.V(i,:));
      end
    end

    % show the 3D plot of the trajectory
    % optional arguments: from=1, to=self.T to splice the range
		% also color the beginning and the end of the track green and red
		% respectively to specify direction of time
    function show(self, varargin)
      from = 1;
      to = self.T;
      R = self.positions;

      if nargin > 1, from = varargin{1}; end
      if nargin > 2, to = varargin{2}; end

			interval = floor((to - from + 1) / 10);
			iBegin = from + interval; % end of the range at the beginning of the track
			iEnd = to - interval;	% start of the range at the end of the track

      figure;
      plot3(R(from:iBegin,1), 	R(from:iBegin,2), 	R(from:iBegin,3), 	'g', ...
						R(iBegin+1:iEnd,1), R(iBegin+1:iEnd,2), R(iBegin+1:iEnd,3), 'b', ...
						R(iEnd+1:to,1), 		R(iEnd+1:to,2), 		R(iEnd+1:to,3), 		'r');
      title(sprintf(['3D particle trajectory from T=%d to T=%d.' ...
										 'Beginning of track (green), end of the track (red).'], ...
											from, to));
    end

		function plotDiffusion(self) 
			figure;
			plot(self.D);
      ylabel('Diffusion Coefficient ([dist]^2/s)')
      xlabel('n (# steps)');
		end

		function plotVelocity(self)
			figure;
			plot(self.V);
			legend('v_x', 'v_y', 'v_z');
      ylabel('Instantaneous Velocities in each direction ([dist]/s)');
      xlabel('n (# steps)');
		end

    % show summary statistics for trajectory between from and to
    function summarize(self, varargin)
      from = 1;
      to = self.T;

      if nargin > 1, from = varargin{1}; end
      if nargin > 2, to = varargin{2}; end

      figure;
      subplot(1,2,1);
      hist(self.steps(from:to-1,:));
      legend('x', 'y', 'z');
      xlabel('Step size (microns)');
      title(sprintf('Step distribution from T=%d to T=%d', from, to));
       
      subplot(1,2,2);
      hist(self.angles(from:to-2));
      xlabel('Angle (radians)');
      title(sprintf('Angle distribution from T=%d to T=%d', from, to));
    end

    % compare: determine the errors from given predictions from
    % true values if known
    function [errs] = compare(self, D, V, useMag)
			if (nargin < 4)
				useMag = false;
			end
			errs(1) = Track.percentError(self.D, D);

			if (useMag)
				errs(2) = Track.percentError(self.normV, V);
			else
				for j=1:3
					errs(j+1) = Track.percentError(self.V(:,j), V(:,j));
				end
			end
    end

    % plus: add another track at the end of this one and return a new instance
    % all tracks implicitly start from a random point in the unit box
    % so to join tracks, we translate the other track to start from the
    % endpoint of this track
    % There is one less point than the sum of points in both tracks, but the
    % number of steps add and the number of angles is one more than the sum in
    % both tracks (there's a new angle at the joining point).
    function self = plus(self, other)
      if self.tau ~= other.tau
        disp(sprintf('Error: two tracks have different timesteps %d vs %d', ...
            self.tau, other.tau));
        return
      end

      last = self.positions(end,:);
      for i=1:other.T-1
        self.positions = [self.positions; last + other.steps(i,:)];
        last = self.positions(end,:);
      end

      % add the angle between the last step of self and first step of other
      self.angles = [self.angles; Track.angle(self.steps(end,:), ...
                                              other.steps(1,:))];

      % join the steps and D,V's for each step
      self.steps = cat(1, self.steps, other.steps);
      self.D = cat(1, self.D, other.D);
      self.V = cat(1, self.V, other.V);

      % now add the rest of the angles from other
      self.angles = cat(1, self.angles, other.angles);

      self.T = self.T + other.T - 1; % other's startpoint is self's endpoint
    end

  end
    
end
