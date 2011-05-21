% Run a suite of test applications on all algorithms.
% May take a while. Also may want to do test in steps to check intermediate
% results (e.g. run bictrain first and then runsolvers) because poor training
% will yield bad final results.

N = 10;			% number of tracks to sample
T = 100;		% length of tracks in steps

diff = 5;		% diffusion step to use in test apps
vel = 5;		% velocity step to use in test apps
dur = 50;		% length of state duration for pulse app

% Suite of test apps
ta = TrackApps(T);
apps(1) = ta.diff();												% DIFF						
apps(2) = ta.diffstep(diff);								% DIFFSTEP
apps(3) = ta.velstep(vel);									% VELSTEP
apps(4) = ta.diffvelstep(diff, vel);				% DIFFVELSTEP

ta.T = 2*T;
apps(5) = ta.pulse(diff, dur); 							% PULSE (need longer tracks)

% Learn cHMM model size
fprintf('Model Selection with BIC\n');
for i=1:5
	fprintf('BIC for app %i: %s\n', i, apps(i).name);
	models(i) = bictrain(N, apps(i));
end

% fHMM parameters: slightly loose max params generally works well
DMAX = diff+1;
VMAX = vel+1;

for i=1:5
	fprintf('Solving app %i: %s\n', i, apps(i).name);
	out(i) = runsolvers(N, apps(i), DMAX, VMAX, models(i).Kstar);
end
