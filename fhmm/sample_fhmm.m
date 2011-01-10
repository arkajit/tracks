% samples nseq length m output sequences from the hmm
function [X, S, V] = sample_fhmm(hmm,nseq,m)

X = cell(nseq,1); 
S = cell(nseq,1);
V = cell(nseq,1);
for i=1:nseq,
    x = zeros(1,m);
    s = zeros(1,m);
    s(1) = sample(hmm.t); 
    [m1, s1] = hmm.s2ind(s(1));    
    x(1) = hmm.means(m1) + randn * hmm.stddevs(s1);
    for j=2:m,
        s(j) = sample(hmm.T(s(j-1),:)); 
        [mj, sj] = hmm.s2ind(s(j));    
        x(j) = hmm.means(mj) + randn * hmm.stddevs(sj);
    end;
    X{i} = x; 
    V{i} = viterbi(x,hmm);
    S{i} = s;
end;


% -----------------
function [s] = sample(p)

p = p/sum(p); % just in case 

r = rand; 
psum = p(1); s=1;
while (psum<r),
    s = s+1;
    psum = psum + p(s); 
end;
