apps(1) = TrackApps.diff();
apps(2) = TrackApps.diffstep(5);
apps(3) = TrackApps.velstep(5);
apps(4) = TrackApps.diffvelstep(5, 5);
apps(5) = TrackApps.pulse(5, 50);

N = 100;
T = 100;

for i=1:4
	fprintf('app %i\n', i);
	out(i) = bictrain(N, T, apps(i));
end

fprintf('pulse\n');
out(5) = bictrain(N, 2*T, apps(5));
