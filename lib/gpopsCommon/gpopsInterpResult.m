function result = gpopsInterpResult(result)

% gpopsInterpResult
% interp the result to a grid of 5*N linear spaced points in each phase

solution = result.solution;

setup = result.setup;

% get number of phases
numphases = size(result.solution.phase,2);

% add 1 to the number of collocation points
currentmesh = setup.mesh.phase;

% preallocate interpsolution
interpphase(numphases).time = [];
interpphase(numphases).state = [];
interpphase(numphases).control = [];
if isfield(solution.phase, 'integral');
    interpphase(numphases).integral = [];
end
interpphase(numphases).costate = [];

% interp solution, control
for phasecount = 1:numphases;
    % get solution for each phase
    phasesol = solution.phase(phasecount);
    
    % get result mesh and interp for each phase
    resultcolpoints = currentmesh(phasecount).colpoints;
    ssol = [result.collocation(phasecount).s; [1 1]];
    
    % number of interp points
    numnodesp = 5*sum(currentmesh(phasecount).colpoints) + 1;
    
    % time on interp mesh for phase
    t0 = phasesol.time(1);
    tf = phasesol.time(end);
    
    % get number of segments
    numseg = size(resultcolpoints,2);
    
    % preallocate time, state, control and costate for each phase
    timeinterp = zeros(numnodesp, 1);
    numstate = size(phasesol.state,2);
    stateinterp = zeros(numnodesp, numstate);
    if isfield(phasesol, 'control');
        if ~isempty(phasesol.control);
            numcontrol = size(phasesol.control,2);
            controlinterp = zeros(numnodesp, numcontrol);
            phasecontrolswitch = true;
        else
            controlinterp = [];
            phasecontrolswitch = false;
        end
    else
        controlinterp = [];
        phasecontrolswitch = false;
    end
    costateinterp = zeros(numnodesp, numstate);
    
    % initalize segment markers
    solstart = 1;
    intstart = 1;
    for segcount = 1:numseg;
        % update end marker
        solend = solstart + resultcolpoints(segcount);
        intend = intstart + 5*resultcolpoints(segcount)-1;
        
        % get solution and interpolation index
        solindex = (solstart:solend)';
        intindex = (intstart:intend)';
        
        % get segment LGR points
        ssolseg = ssol(solindex,2);
        ssolseg(end) = 1;
        sinterpseg = linspace(-1, 1, 5*resultcolpoints(segcount)+1)';
        sinterpseg = sinterpseg(1:5*resultcolpoints(segcount));
        
        % interp time using lagrange polynomial
        segt0 = phasesol.time(solstart);
        segtf = phasesol.time(solend);
        timeinterp(intindex) = (sinterpseg + 1).*(segtf - segt0)./2 + segt0;
        
        % interp state using lagrange polynomial
        segstate = phasesol.state(solindex,:);
        stateinterp(intindex,:) = gpopsLagrangeInterp(ssolseg, segstate, sinterpseg);
        
        % interp control using cubic
        if phasecontrolswitch;
            segcontrol = phasesol.control(solindex,:);
            
            %interp using cubic
            controlinterp(intindex,:) = interp1(ssolseg, segcontrol, sinterpseg, 'pchip');
        end
        
        % update start marker
        solstart = solend;
        intstart = intend+1;
    end
    % Get full state (state at radau points including end point at 1
    timeinterp(intstart) = tf;
    stateinterp(intstart,:) = phasesol.state(end,:);
    if phasecontrolswitch;
        controlinterp(intstart,:) = phasesol.control(end,:);
    end
    costateinterp(intstart,:) = phasesol.costate(end,:);
    
    % save interp state and control for each phase
    interpphase(phasecount).time = timeinterp;
    interpphase(phasecount).state = stateinterp;
    interpphase(phasecount).control = controlinterp;
    if isfield(solution.phase, 'integral');
        interpphase(phasecount).integral = solution.phase(phasecount).integral;
    end
    interpphase(phasecount).costate = costateinterp;
end

interpsolution.phase = interpphase;

if isfield(solution, 'parameter');
    interpsolution.parameter = solution.parameter;
end

result.interpsolution = interpsolution;