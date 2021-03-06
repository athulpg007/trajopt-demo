function [continput, endpinput, tp0, tpf] = gpopsContEndpInputRPMI(Z, probinfo)

% gpopsContEndpInputRPMI
% this function gets the inputs for the continuous and endpoint functions
% from the NLP variable vector
% tp0 and tpf contain the initial and final time for each phase

% get OCP info
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numintegral = probinfo.numintegral;
numparameter = probinfo.numparameter;

% get number of nodes
numnodes = probinfo.numnodes;

% preallocate contphase and endpphase
contphase(numphase).state = [];
endpphase(numphase).initialstate = [];
endpphase(numphase).finalstate = [];
if sum(numcontrol,2) ~= 0;
    contphase(numphase).control = [];
end
contphase(numphase).time = [];
endpphase(numphase).initialtime = [];
endpphase(numphase).finaltime = [];
if numparameter ~= 0;
    contphase(numphase).parameter = [];
end
if sum(numintegral,2) ~= 0;
    endpphase(numphase).integral = [];
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
    numintegralp = numintegral(phasecount);
    numnodesp = numnodes(phasecount);
    
    % get nlp map for phase
    phasenlpmap = probinfo.nlpcontmap(phasecount);
    
    % get OCP state for phase
    state = Z(phasenlpmap.statemap(1,1):phasenlpmap.statemap(2,numstatep));
    state = reshape(state,numnodesp+1,numstatep);
    contphase(phasecount).state = state(1:numnodesp,:);
    endpphase(phasecount).initialstate = state(1,:);
    endpphase(phasecount).finalstate   = state(numnodesp+1,:);
    
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
    endpphase(phasecount).initialtime = t0;
    endpphase(phasecount).finaltime   = tf;
    
    % get OCP integral for phase
    if numintegralp ~= 0;
        endpphase(phasecount).integral = Z(phasenlpmap.integralvarmap)';
    end
    
    % size parameters for continput in each phase
    if numparameter ~= 0;
        contphase(phasecount).parameter = ones(numnodesp,1)*parameter;
    end
end

% add variables for all phases to continput and enpinput
continput.phase = contphase;
endpinput.phase = endpphase;

% get endpinput for parameter guess
if numparameter ~= 0;
    endpinput.parameter = parameter;
end

% add auxdata to continput and endpinput
if probinfo.auxflag;
    continput.auxdata = probinfo.auxdata;
    endpinput.auxdata = probinfo.auxdata;
end