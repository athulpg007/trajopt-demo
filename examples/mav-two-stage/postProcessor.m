solution = output.result.solution;        % load solution structure 
tf       = solution.phase(2).time(end);   % extract final time tf
trem     = 569 - tf;                      % stage 2 remaining burn time
frem     = trem*mdotSecond;               % stage 2 fuel mass remaining

mf0      = solution.phase(2).state(end,7); % extract final mass
mff      = mf0 - frem ;                    % final mass when 2nd stage is used completely

dV       = log(mf0/mff)*360*9.80665

 