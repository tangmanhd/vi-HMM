function [Delta, xtilde] = calculateDeltaHMM(mysequence, initialprobs, transitionmatrix, emissionmatrix)
% mysequence: vector of length L
%initialprobs: vector of length M, M: number of possible states
% transitionmatrix: matrix of size M by M
% emissionmatrixful: matrix of size M by N, N: number of possible observations
L = length(mysequence);
M = length(initialprobs);
Delta = NaN(M, L);
xtilde = NaN(1, L);

%Delta(:, 1) = initialprobs(:, 1) .* emissionmatrix(:,1);
Delta(:, 1) = initialprobs(:, 1) .* exp(emissionmatrix(:,1));
Delta(:, 1) = Delta(:, 1) / sum(Delta(:, 1));
for t = 2 : L
%     usually do this in log scale
    Delta(:, t) = log(max(repmat(Delta(:, t - 1), [1, M]) .* transitionmatrix)' )+ emissionmatrix(:,t);
    den=logsumexp(Delta(:, t));
    Delta(:, t)=exp(Delta(:, t)-den);
   % Delta(:, t) = max(repmat(Delta(:, t - 1), [1, M]) .* transitionmatrix)' .* emissionmatrix(:,t);
   % Delta(:, t) = Delta(:, t) / sum(Delta(:, t));
end
xtilde(L) = find(Delta(:, L) == max(Delta(:, L)));
for t = (L - 1) : -1 : 1
    temp = Delta(:, t) .* transitionmatrix(:, xtilde(t + 1));
    xtilde(t) = find(temp == max(temp),1);
end
end

