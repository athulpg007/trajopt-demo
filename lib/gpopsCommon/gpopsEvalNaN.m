function [contdependnan, endpdependnan] = gpopsEvalNaN(setup, probinfo)

% gpopsEvalNaN
% this function is used to find the dependencies of the optimal control
% problem, the optimal control problem is evaluated on the initial guess
% where a single variable is given the value 'NaN', the resulting optimal 
% control problem output is checked for values of 'NaN', the functions with
% the values of 'NaN' are then marked as dependent on the variable that was
% given the value 'NaN', this process is repeated for all variables in the
% optimal control prolblem

% get OCP sizes
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numpath = probinfo.numpath;
numintegral = probinfo.numintegral;
numparameter = probinfo.numparameter;

numeventgroup = probinfo.numeventgroup;
numevent = probinfo.numevent;

% get inputs for continuous and endpoint functions
% preallocate contphase and endpphase
contphase(numphase).state = [];
endpphase(numphase).initialstate = [];
endpphase(numphase).finalstate = [];
if sum(numcontrol,2) ~= 0;
    controlswitch = true;
    contphase(numphase).control = [];
else
    controlswitch = false;
end
contphase(numphase).time = [];
endpphase(numphase).initialtime = [];
endpphase(numphase).finaltime = [];
if numparameter ~= 0;
    contphase(numphase).parameter = [];
end
if sum(numintegral,2) ~= 0;
    integralswitch = true;
    endpphase(numphase).integral = [];
else
    integralswitch = false;
end

guess = setup.guess;

% get input values for each phase
for phasecount = 1:numphase;
    % get guess for each phase
    solphase = guess.phase(phasecount);
    
    % get input for state guess
    contphase(phasecount).state = solphase.state;
    endpphase(phasecount).initialstate = solphase.state(1,:);
    endpphase(phasecount).finalstate = solphase.state(end,:);
    
    % get input for control guess
    if controlswitch;
        contphase(phasecount).control = solphase.control;
    end
    
    % get input for time guess
    contphase(phasecount).time = solphase.time;
    endpphase(phasecount).initialtime = solphase.time(1);
    endpphase(phasecount).finaltime = solphase.time(end);
    
    % get contphase for parameter guess
    if numparameter ~= 0;
        contphase(phasecount).parameter = ones(length(solphase.time),1)*guess.parameter;
    end
    
    % get endpphase for integral guess
    if integralswitch;
        endpphase(phasecount).integral = solphase.integral;
    end
end

% adding guess for all phases to continput and endpinput
continput.phase = contphase;
endpinput.phase = endpphase;

% get endpinput for parameter guess
if numparameter ~= 0;
    endpinput.parameter = guess.parameter;
end

% add auxdata to continput and endpinput
if isfield(setup,'auxdata');
    continput.auxdata = setup.auxdata;
    endpinput.auxdata = setup.auxdata;
end

% find continuous function dependencies -----------------------------------
% get continuous function derivative map
contmap1 = probinfo.derivativemap.contmap1;

% get nvc for first derivatives
contnvc = probinfo.derivativemap.contnvc1;
maxnvc = max(contnvc);

% preallocate dependency output for continuous function
dynamicsnnz = probinfo.derivativemap.dynamicsnnz1;
contdependnan(numphase).dynamicsgrd = [];
if sum(numpath,2) ~= 0;
    pathnnz = probinfo.derivativemap.pathnnz1;
    contdependnan(numphase).pathgrd = [];
else
    pathnnz = zeros(1,numphase);
end
if integralswitch;
    integrandnnz = probinfo.derivativemap.integrandnnz1;
    contdependnan(numphase).integrandgrd = [];
else
    integrandnnz = zeros(1,numphase);
end
for phasecount = 1:numphase;
    NN = size(guess.phase(phasecount).time,1);
    if dynamicsnnz(phasecount) ~= 0;
        contdependnan(phasecount).dynamicsgrd = zeros(NN,dynamicsnnz(phasecount));
    end
    if pathnnz(phasecount) ~= 0;
        contdependnan(phasecount).pathgrd = zeros(NN,pathnnz(phasecount));
    end
    if integrandnnz(phasecount) ~= 0;
        contdependnan(phasecount).integrandgrd = zeros(NN,integrandnnz(phasecount));
    end
end

% find continuous function dependencies of all phases simultaneously
for nvccount = 1:maxnvc;
    % initiate pertinput as the unperturbed input
    pertinput = continput;
    for phasecount = 1:numphase;
        numstatep = numstate(phasecount);
        numcontrolp = numcontrol(phasecount);
        % test if the variable count is less then the number of
        % variables in the phase
        if nvccount <= contnvc(phasecount);
            varnum = contmap1(phasecount).contvarmap1(nvccount);
            if varnum <= numstatep;
                % state
                refmark = varnum;
                pertinput.phase(phasecount).state(:,refmark) = NaN;
            elseif varnum <= numstatep+numcontrolp;
                % control
                refmark = varnum-numstatep;
                pertinput.phase(phasecount).control(:,refmark) = NaN;
            elseif varnum == numstatep+numcontrolp+1;
                % time
                pertinput.phase(phasecount).time(:) = NaN;
            elseif varnum <= numstatep+numcontrolp+1+numparameter;
                % parameter
                refmark = varnum-numstatep-numcontrolp-1;
                pertinput.phase(phasecount).parameter(:,refmark) = NaN;
            end
        end
    end
    
    % evaluate function on perturbed input
    pertoutput = feval(probinfo.contfunction, pertinput);
    
    % find the dependencies in each phase
    for phasecount = 1:numphase;
        numstatep = numstate(phasecount);
        numpathp = numpath(phasecount);
        numintegralp = numintegral(phasecount);
        % test if the perturbed variable count is less then the number of
        % perturbed variables in the phase
        if nvccount <= contnvc(phasecount);
            % find dependencies of dynamic constraints
            for dynamicscount = 1:numstatep;
                refmark = contmap1(phasecount).dynamicsmap1(dynamicscount,nvccount);
                if refmark ~= 0;
                    contdependnan(phasecount).dynamicsgrd(:,refmark) = isnan(pertoutput(phasecount).dynamics(:,dynamicscount));
                end
            end
            % find dependencies of path constraints
            for pathcount = 1:numpathp;
                refmark = contmap1(phasecount).pathmap1(pathcount,nvccount);
                if refmark ~= 0;
                    contdependnan(phasecount).pathgrd(:,refmark) = isnan(pertoutput(phasecount).path(:,pathcount));
                end
            end
            % find dependencies of intergrand constraints
            for intergralcount = 1:numintegralp;
                refmark = contmap1(phasecount).integrandmap1(intergralcount,nvccount);
                if refmark ~= 0;
                    contdependnan(phasecount).integrandgrd(:,refmark) = isnan(pertoutput(phasecount).integrand(:,intergralcount));
                end
            end
        end
    end
end

% find endpoint dependencies ----------------------------------------------

% get endpoint derivative map
endpvarloc = probinfo.endpvarloc;
endpnvc1 = probinfo.derivativemap.endpnvc1;
objnnz1 = probinfo.derivativemap.objnnz1;
endpvarmap1 = probinfo.derivativemap.endpvarmap1;
endpobjmap1 = probinfo.derivativemap.endpobjmap1;

% preallocate endpoint dependency output
objectivegrd = zeros(1,objnnz1);
if numeventgroup ~= 0;
    eventnnz1 = probinfo.derivativemap.eventnnz1;
    endpeventmap1 = probinfo.derivativemap.endpeventmap1;
    eventgroup(numeventgroup).eventgrd = [];
    for eventgroupcount = 1:numeventgroup;
        eventgroup(eventgroupcount).eventgrd = zeros(1,eventnnz1(eventgroupcount));
    end
end

for nvccount = 1:endpnvc1
    % initiate pertinput as the unperturbed input
    pertinput = endpinput;
    
    % get variable phase, and variable number in phase
    varnum = endpvarmap1(nvccount);
    varphase = endpvarloc(1,varnum);
    phasevarnum = endpvarloc(2,varnum);
    
    if varphase ~= 0;
        if phasevarnum <= numstate(varphase);
            % initial state
            refmark = phasevarnum;
            pertinput.phase(varphase).initialstate(refmark) = NaN;
        elseif phasevarnum <= 2*numstate(varphase);
            % final state
            refmark = phasevarnum-numstate(varphase);
            pertinput.phase(varphase).finalstate(refmark) = NaN;
        elseif phasevarnum == 2*numstate(varphase)+1;
            % initial time
            pertinput.phase(varphase).initialtime = NaN;
        elseif phasevarnum == 2*numstate(varphase)+2;
            % final time
            pertinput.phase(varphase).finaltime = NaN;
        elseif phasevarnum <= 2*numstate(varphase)+2+numintegral(varphase);
            % integral
            refmark = phasevarnum-2*numstate(varphase)-2;
            pertinput.phase(varphase).integral(refmark) = NaN;
        end
    else
        % parameter
        refmark = phasevarnum;
        pertinput.parameter(refmark) = NaN;
    end
    
    % evaluate function on perturbed input
    pertoutput = feval(probinfo.endpfunction, pertinput);
    
    % find dependencies of the objective
    refmark = endpobjmap1(nvccount);
    if refmark ~= 0;
        objectivegrd(refmark) = isnan(pertoutput.objective);
    end
    
    % find dependencies of event constraints
    for eventgroupcount = 1:numeventgroup;
        for eventcount = 1:numevent(eventgroupcount);
            refmark = endpeventmap1(eventgroupcount).first(eventcount,nvccount);
            if refmark ~= 0;
                eventgroup(eventgroupcount).eventgrd(refmark) = isnan(pertoutput.eventgroup(eventgroupcount).event(eventcount));
            end
        end
    end
end

endpdependnan.objectivegrd = objectivegrd;
if numeventgroup ~= 0;
    endpdependnan.eventgroup = eventgroup;
end
