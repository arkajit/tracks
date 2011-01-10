function [scores, errs, etypes] = fact()
  means = [0; 5; 10];
  stddevs = [1; 2; 3];
  t = [1; 0; 0];
  odds=1:100;
  scores = zeros(1,length(odds));
  errs = cell(1,length(odds));
  etypes = zeros(9);

  for o=odds
    f = fact_hmm(means, stddevs, t, t, o, o);
    [X, S, V] = sample_fhmm(f,10,100);
    [scores(o), errs{o}] = score(f,S,V); 
     for e=1:size(errs{o},1)
       i = errs{o}(e,:);
       etypes(i(1),i(2)) = etypes(i(1),i(2))+1;
     end
  end

  plot(scores/2000); % 2000 = 10*(2*100) number of states to predict
  xlabel('Odds against state switching');
  ylabel('% of states that differ between Viterbi and true paths');
end

function [total, errors] = score(f,S,V)
  n = length(S);
  scores = zeros(1,n);
  errors = [];
  for i=1:n
    diff = V{i}-S{i};
    errs = find(diff); % indices where V and S mismatch
    scores(i) = scores(i) + length(errs);
  
    %penalize extra for states mismatching in both mean and stddev
    for e=errs
      errors = [errors; S{i}(e), V{i}(e)]; %#ok<AGROW>
      [i1,j1] = f.s2ind(errors(end,1));
      [i2,j2] = f.s2ind(errors(end,2));
      if i1 ~= i2 && j1 ~= j2, scores(i) = scores(i) + 1; end;
    end

    %pos_mod = find(mod(diff(errs),3));
    %cores(i) = scores(i) + length(pos_mod); 
      
  end;
  total = sum(scores);
end
