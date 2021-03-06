function yinterp = gpopsLagrangeInterp(s, y, sinterp)

% gpopsLagrangeInterp
% This function interpolates using Lagrange polynomial to find yinterp
% at the values sinterp
% s and sinterp must be on the same domain
% s and y must me the same length

% get number of states and interpolation points
numpoints = size(s,1);
numstate = size(y,2);
numinterp = size(sinterp,1);

% initialize yinterp
yinterp = zeros(numinterp, numstate);
for k = 1:numpoints
    % initialize weights
    w = ones(numinterp,1);
    % for j not equal to k
    for j = [1:k-1 k+1:numpoints]
        w = (sinterp-s(j))./(s(k)-s(j)).*w;
    end
    for statecount = 1:numstate;
        yinterp(:,statecount) = yinterp(:,statecount) + w*y(k,statecount);
    end
end