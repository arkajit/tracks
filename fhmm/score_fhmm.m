function [log_p] = score_fhmm(hmm, s, x)
  log_p = hmm.log_t(s(1)) + hmm.log_E(s(1), x(1));
  for j=2:length(x)
    log_p = log_p + hmm.log_T(s(j-1), s(j)) + hmm.log_E(s(j), x(j));
  end;
