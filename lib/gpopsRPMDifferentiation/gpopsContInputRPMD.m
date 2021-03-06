function [continput, tp0, tpf] = gpopsContInputRPMD(Z, probinfo)

% gpopsContInputRPMD
% this function gets the input for the continuous function
% from the NLP variable vector
% tp0 and tpf contain the initial and final time for each phase

% get OCP info
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numparameter = probinfo.numparameter;

% get number of nodes
numnodes = probinfo.numnodes;

% preallocate contphase
contphase(numphase).state = [];
if sum(numcontrol,2) ~= 0;
    contphase(numphase).control = [];
end
contphase(numphase).time = [];
if numparameter ~= 0;
    contphase(numphase).parameter = [];
end

% get parameters
if numparameter ~= 0;
    parameter = Z(probinfo.nlpparametermap)';
end

% preallocate tp0 and tpf
tp0 = zeros(1,numphase)*Z(1,1);
tpf = zeros(1,numphase)*Z(1,1);

% get variables for phases
for phasecount = 1:numphase;
    % get OCP info for each phase
    numstatep = numstate(phasecount);
    numcontrolp = numcontrol(phasecount);
    numnodesp = numnodes(phasecount);
    
    % get nlp map for phase
    phasenlpmap = probinfo.nlpcontmap(phasecount);
    
    % get OCP state for phase
    state = Z(phasenlpmap.statemap(1,1):phasenlpmap.statemap(2,numstatep));
    state = reshape(state,numnodesp+1,numstatep);
    contphase(phasecount).state = state(1:numnodesp,:);
    
    % get OCP control for phase
    if numcontrolp ~= 0;
        control = Z(phasenlpmap.controlmap(1,1):phasenlpmap.controlmap(2,numcontrolp));
        contphase(phasecount).control = reshape(control,numnodesp,numcontrolp);
    end
    
    % get OCP time for phase
    s = probinfo.collocation(phasecount).s(:,1);
    t0 = Z(phasenlpmap.timemap(1));
    tf = Z(phasenlpmap.timemap(2));
    time = (s + 1).*(tf - t0)./2 + t0;
    tp0(phasecount) = t0;
    tpf(phasecount) = tf;
    contphase(phasecount).time = time;
    
    % size parameters for continput in each phase
    if numparameter ~= 0;
        contphase(phasecount).parameter = ones(numnodesp,1)*parameter;
    end
end

% add variables for all phases to continput and enpinput
continput.phase = contphase;

% add auxdata to continput and endpinput
if probinfo.auxflag;
    continput.auxdata = probinfo.auxdata;
end