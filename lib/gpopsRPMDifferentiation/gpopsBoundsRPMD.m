function [ZL, ZU, FL, FU, probinfo] = gpopsBoundsRPMD(setup, probinfo)

% gpopsBoundsRPMD
% this function defines the variable and constraint bounds for the
% nonlinear program resulting from the Radau Pseudospectral Method
% it also creates a map of optimal control problem variables and
% constraints to the nonlinear program variables and constraints
%
% ZL and ZU represent the lower and upper bounds of the NLP variables
% FL and FU represent the lower and upper bounds of the NLP constraints

% NLP variable order
% [states*(nodes+1); controls*nodes; t0; tf, Q] for each phase
% [stack all phases(1;...;numphase); parameters]

% NLP constraint order
% [defect*nodes; path*nodes; integral; duration] for each phase
% [stack all phases(1;...;numphase); stack all events(1;...;numeventgroup]

% get OCP info
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numpath = probinfo.numpath;
numintegral = probinfo.numintegral;
phaseduration = probinfo.phaseduration;
numparameter = probinfo.numparameter;
numeventgroup = probinfo.numeventgroup;
numevent = probinfo.numevent;

% get number of nodes
numnodes = probinfo.numnodes;

% NLP variable and constraint sizes
nlpnumvar = numstate*(numnodes+1)' + numcontrol*numnodes' + 2*numphase + sum(numintegral,2) + numparameter;
nlpnumcon = (numstate+numpath)*numnodes' + sum(numintegral,2) + numphase + sum(numevent,2);

% preallocate NLP bounds
ZL = zeros(nlpnumvar,1);
ZU = zeros(nlpnumvar,1);
FL = zeros(nlpnumcon,1);
FU = zeros(nlpnumcon,1);

% preallocate NLP map
nlpcontmap(numphase).statemap = [];
if sum(numcontrol,2) ~= 0;
    nlpcontmap(numphase).controlmap = [];
end
nlpcontmap(numphase).timemap = [];
nlpcontmap(numphase).defectmap = [];
if sum(numpath,2) ~= 0;
    nlpcontmap(numphase).pathmap = [];
end
if sum(numintegral,2) ~= 0;
    nlpcontmap(numphase).integralvarmap = [];
    nlpcontmap(numphase).integrandmap = [];
end
nlpcontmap(numphase).durationmap = [];
nlpendpvarmap = zeros(1,length(probinfo.endpvarloc));
if numeventgroup ~= 0;
    nlpeventmap = zeros(2,numeventgroup);
end

% define bounds and NLP map
Zasgmark = 0;
Fasgmark = 0;
endpvarmark = 0;
for phasecount = 1:numphase;
    % OCP info for phase
    numstatep = numstate(phasecount);
    numcontrolp = numcontrol(phasecount);
    numpathp = numpath(phasecount);
    numintegralp = numintegral(phasecount);
    numnodesp = numnodes(phasecount);
    phasedurationp = phaseduration(phasecount);
    
    % get phase bounds
    phasebounds = setup.bounds.phase(phasecount);
    
    % preallocate statemap and dynamicsmap for each phase
    statemap = zeros(2,numstatep);
    defectmap = zeros(2,numstatep);
    
    % define bounds for state variables and defect bounds
    for statecount = 1:numstatep;
        % get Z and F index start and end points
        Zindexs = Zasgmark + 1;
        Zindexe = Zasgmark + (numnodesp+1);
        Findexs = Fasgmark + 1;
        Findexe = Fasgmark + numnodesp;
        
        % set state and defect map
        statemap(:,statecount) = [Zindexs; Zindexe];
        defectmap(:,statecount) = [Findexs; Findexe];
        
        % initial state location
        nlpendpvarmap(endpvarmark+statecount) = Zindexs;
        
        % final state location
        nlpendpvarmap(endpvarmark+numstatep+statecount) = Zindexe;
        
        % FL and FU already zeros for these locations
        % set ZL and ZU
        ZL(Zindexs) = phasebounds.initialstate.lower(statecount);
        ZU(Zindexs) = phasebounds.initialstate.upper(statecount);
        Zindex = Zindexs+1:Zindexe-1;
        ZL(Zindex) = phasebounds.state.lower(statecount);
        ZU(Zindex) = phasebounds.state.upper(statecount);
        ZL(Zindexe) = phasebounds.finalstate.lower(statecount);
        ZU(Zindexe) = phasebounds.finalstate.upper(statecount);
        
        % change marker value
        Zasgmark = Zindexe;
        Fasgmark = Findexe;
    end
    % change endp variable marker
    endpvarmark = endpvarmark + 2*numstatep;
    
    % define state map and defect map for phase
    nlpcontmap(phasecount).statemap = statemap;
    nlpcontmap(phasecount).defectmap = defectmap;
    
    % check if phase contains control variables
    if numcontrolp ~= 0;
        % preallocate controlmap
        controlmap = zeros(2,numcontrolp);
        
        % define bounds for control variables
        for controlcount = 1:numcontrolp;
            % get Z index start and end points
            Zindexs = Zasgmark + 1;
            Zindexe = Zasgmark + numnodesp;
            
            % set control map
            controlmap(:,controlcount) = [Zindexs; Zindexe];
            
            % set ZL and ZU
            Zindex = Zindexs:Zindexe;
            ZL(Zindex) = phasebounds.control.lower(controlcount);
            ZU(Zindex) = phasebounds.control.upper(controlcount);
            
            % change marker value
            Zasgmark = Zindexe;
        end
        % define control map for phase
        nlpcontmap(phasecount).controlmap = controlmap;
    end
    % check if phase contains path constraints
    if numpathp ~= 0;
        % preallocate pathmap
        pathmap = zeros(2,numpathp);
        
        % define bounds for path constraints
        for pathcount = 1:numpathp;
            % get Z and F index start and end points
            Findexs = Fasgmark + 1;
            Findexe = Fasgmark + numnodesp;
            
            % set path map
            pathmap(:,pathcount) = [Findexs; Findexe];
            
            % set FL and FU
            Findex = Findexs:Findexe;
            FL(Findex) = phasebounds.path.lower(pathcount);
            FU(Findex) = phasebounds.path.upper(pathcount);
            
            % change marker value
            Fasgmark = Findexe;
        end
        % define path map for phase
        nlpcontmap(phasecount).pathmap = pathmap;
    end
    % define bounds for initial time variable
    Zasgmark = Zasgmark + 1;
    ZL(Zasgmark) = phasebounds.initialtime.lower;
    ZU(Zasgmark) = phasebounds.initialtime.upper;
    
    % initial time location
    nlpendpvarmap(endpvarmark+1) = Zasgmark;
    nlpcontmap(phasecount).timemap(1,1) = Zasgmark;
    
    % define bounds for final time variable
    Zasgmark = Zasgmark + 1;
    ZL(Zasgmark) = phasebounds.finaltime.lower;
    ZU(Zasgmark) = phasebounds.finaltime.upper;
    
    % final time location
    nlpendpvarmap(endpvarmark+2) = Zasgmark;
    nlpcontmap(phasecount).timemap(2,1) = Zasgmark;
    
    % change endp variable marker
    endpvarmark = endpvarmark + 2;
    
    % check if phase contains integral constraints
    if numintegralp ~= 0;
        % preallocate integral map
        integralvarmap = zeros(1,numintegralp);
        integrandmap = zeros(1,numintegralp);
        for integralcount = 1:numintegralp;
            % get Z and F location
            Zasgmark = Zasgmark + 1;
            Fasgmark = Fasgmark + 1;
            
            % set integral map
            integralvarmap(integralcount) = Zasgmark;
            integrandmap(integralcount) = Fasgmark;
            
            % integral variable location
            nlpendpvarmap(endpvarmark+integralcount) = Zasgmark;
            
            % set ZL and ZU (FL and FU already zero)
            ZL(Zasgmark) = phasebounds.integral.lower(integralcount);
            ZU(Zasgmark) = phasebounds.integral.upper(integralcount);
        end
        % define integral map for phase
        nlpcontmap(phasecount).integralvarmap = integralvarmap;
        nlpcontmap(phasecount).integrandmap = integrandmap;
        
        % change endp variable marker
        endpvarmark = endpvarmark + numintegralp;
    end
    
    % check if duration is user defined in phase
    Fasgmark = Fasgmark + 1;
    if phasedurationp;
        % set user defined duration bounds dU >= tf - t0 >= dL
        FL(Fasgmark) = phasebounds.duration.lower;
        FU(Fasgmark) = phasebounds.duration.upper;
    else
        % duration bounds so that 2*(tfmax-t0min) >= tf - t0 >= 0
        FL(Fasgmark) = 0;
        FU(Fasgmark) = 2*(phasebounds.finaltime.upper - phasebounds.initialtime.lower);
    end
    nlpcontmap(phasecount).durationmap = Fasgmark;
end

% check if problem has events
if numeventgroup ~= 0;
    % get event bounds
    eventbounds = setup.bounds.eventgroup;
    for eventgroupcount = 1:numeventgroup;
        % get F index points
        Findexs = Fasgmark + 1;
        Findexe = Fasgmark + numevent(eventgroupcount);
        Findex = Findexs:Findexe;
        
        % define event map
        nlpeventmap(1,eventgroupcount) = Findexs;
        nlpeventmap(2,eventgroupcount) = Findexe;
        
        % define event bounds
        FL(Findex) = eventbounds(eventgroupcount).lower;
        FU(Findex) = eventbounds(eventgroupcount).upper;
        
        % change marker value
        Fasgmark = Findexe;
    end
end

% check is problem has parameters
if numparameter ~= 0;
    % get Z index points
    Zindexs = Zasgmark + 1;
    Zindexe = Zasgmark + numparameter;
    Zindex = Zindexs:Zindexe;
    
    % define parameter map
    nlpparametermap = Zindex;
    
    % define parameter location in endpvarmap
    nlpendpvarmap((1:numparameter)+endpvarmark) = Zindex;
    
    % define parameter bounds
    ZL(Zindex) = setup.bounds.parameter.lower;
    ZU(Zindex) = setup.bounds.parameter.upper;
end

% NLP sizes to probinfo
probinfo.nlpnumvar = nlpnumvar;
probinfo.nlpnumcon = nlpnumcon;

% add NLP map to probinfo
% preallocate NLP map
probinfo.nlpcontmap = nlpcontmap;
if numparameter ~= 0;
    probinfo.nlpparametermap = nlpparametermap;
end
probinfo.nlpendpvarmap = nlpendpvarmap;
if numeventgroup ~= 0;
    probinfo.nlpeventmap = nlpeventmap;
end