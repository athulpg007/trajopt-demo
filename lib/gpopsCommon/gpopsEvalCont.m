function contoutput = gpopsEvalCont(solution, setup)

% gpopsEvalCont
% this function evaluates the continuous function on the solution

% get number of phases
numphase = length(solution.phase);

% get number of parameters
if isfield(solution,'parameter');
    numparameter = length(solution.parameter);
else
    numparameter = 0;
end

% preallocate contphase and endpphase
contphase(numphase).state = [];
if isfield(solution.phase, 'control');
    controlswitch = true;
    contphase(numphase).control = [];
else
    controlswitch = false;
end
contphase(numphase).time = [];
if numparameter ~= 0;
    contphase(numphase).parameter = [];
end

% get input values for each phase
for phasecount = 1:numphase;
    % get solution for each phase
    solphase = solution.phase(phasecount);
    
    % get input for state solution
    contphase(phasecount).state = solphase.state;
    
    % get input for control solution
    if controlswitch;
        contphase(phasecount).control = solphase.control;
    end
    
    % get input for time solution
    contphase(phasecount).time = solphase.time;
    
    % get contphase for parameter solution
    if numparameter ~= 0;
        contphase(phasecount).parameter = ones(length(solphase.time),1)*solution.parameter;
    end
end

% adding solution for all phases to continput
continput.phase = contphase;

% add auxdata to continput and endpinput
if isfield(setup,'auxdata');
    continput.auxdata = setup.auxdata;
end

% evaluate OCP continuous function
contoutput = feval(setup.functions.continuous, continput);