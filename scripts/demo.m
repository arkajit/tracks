% Compute BIC scores for all test applications and plot results

N = 100;		% number of tracks to sample
T = 100;		% length of tracks in steps

ta = TrackApps(T);
apps(1) = ta.diff();
apps(2) = ta.diffstep(5);
apps(3) = ta.velstep(5);
apps(4) = ta.diffvelstep(5, 5);
apps(5) = TrackApps(2*T).pulse(5, 50); % need longer tracks for good results

for i=1:5
	fprintf('app %i\n', i);
	out(i) = bictrain(N, apps(i));
end
