function [contgrdsamples, endpgrdsamples] = gpopsRandomGrd(setup, probinfo, numsamples)

% gpopsRandomGrd
% this function computes the first derivatives of the optimal control
% problem using the selected finite difference method at random points
% within the variable bounds, the number of random points is determined by
% the numsamples input

% set rand state so that random numbers are repeatable
if exist('rng') > 0;
    rng('default');
else
    if exist('randstream') > 0;
        s = RandStream.create('mt19937ar','seed',5489);
        RandStream.setDefaultStream(s);
    else
        rand('state',0)
    end
end

% generate random sample points -------------------------------------------
% get OCP sizes
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numintegral = probinfo.numintegral;
numparameters = probinfo.numparameter;
numeventgroup = probinfo.numeventgroup;
sampleones = ones(numsamples,1);

% preallocate cphaserand
% cphaserand stores the random points for the continuous function
cphaserand(numphase).state = [];
if sum(numcontrol,2) ~= 0;
    cphaserand(numphase).control = [];
end
cphaserand(numphase).time = [];

% preallocate ephaserand
% ephaserand stores the random points for the endpoiont function
ephaserand(numphase).initialstate = [];
ephaserand(numphase).finalstate = [];
ephaserand(numphase).initialtime = [];
ephaserand(numphase).finaltime = [];
if sum(numintegral,2) ~= 0;
    ephaserand(numphase).integral = [];
end

% get parameter points
if numparameters ~= 0;
    cphaserand(numphase).parameter = [];
    paralower = sampleones*setup.bounds.parameter.lower;
    paraupper = sampleones*setup.bounds.parameter.upper;
    pararand = paralower+(paraupper-paralower).*rand(numsamples,numparameters);
    for phasecount = 1:numphase;
        cphaserand(phasecount).parameter = pararand;
    end
else
    pararand = [];
end

phasebounds = setup.bounds.phase;
% get random cont input values
for phasecount = 1:numphase;
    % get OCP info for each phase
    numstatep = numstate(phasecount);
    numcontrolp = numcontrol(phasecount);
    numintegralp = numintegral(phasecount);
    
    % get state bounds
    statelower = sampleones*phasebounds(phasecount).state.lower;
    stateupper = sampleones*phasebounds(phasecount).state.upper;
    istatelower = sampleones*phasebounds(phasecount).initialstate.lower;
    istateupper = sampleones*phasebounds(phasecount).initialstate.upper;
    fstatelower = sampleones*phasebounds(phasecount).finalstate.lower;
    fstateupper = sampleones*phasebounds(phasecount).finalstate.upper;
    
    % random state points
    cphaserand(phasecount).state = statelower+(stateupper-statelower).*rand(numsamples,numstatep);
    ephaserand(phasecount).initialstate = istatelower+(istateupper-istatelower).*rand(numsamples,numstatep);
    ephaserand(phasecount).finalstate = fstatelower+(fstateupper-fstatelower).*rand(numsamples,numstatep);
    
    if numcontrolp ~= 0;
        % get control bounds
        controllower = sampleones*phasebounds(phasecount).control.lower;
        controlupper = sampleones*phasebounds(phasecount).control.upper;
        
        % random control points
        cphaserand(phasecount).control = controllower+(controlupper-controllower).*rand(numsamples,numcontrolp);
    end
    
    % get time bounds
    itimelower = sampleones*phasebounds(phasecount).initialtime.lower;
    itimeupper = sampleones*phasebounds(phasecount).initialtime.upper;
    ftimelower = sampleones*phasebounds(phasecount).finaltime.lower;
    ftimeupper = sampleones*phasebounds(phasecount).finaltime.upper;
    
    % random time points
    cphaserand(phasecount).time = itimelower+(ftimeupper-itimelower).*sort(rand(numsamples,1));
    ephaserand(phasecount).initialtime = itimelower+(itimeupper-itimelower).*rand(numsamples,1);
    ephaserand(phasecount).finaltime = ftimelower+(ftimeupper-ftimelower).*rand(numsamples,1);
    
    if numintegralp ~= 0;
        integrallower = sampleones*phasebounds(phasecount).integral.lower;
        integralupper = sampleones*phasebounds(phasecount).integral.upper;
        
        ephaserand(phasecount).integral = integrallower+(integralupper-integrallower).*rand(numsamples,numintegralp);
    end
end
% -------------------------------------------------------------------------

% evaluate cont function at random points ---------------------------------
% get continuous function input, evaluate all sample points at once
continput.phase = cphaserand;
% add auxiliary data to continput structure
if probinfo.auxflag;
    continput.auxdata = probinfo.auxdata;
end

% evaluate continuous function derivative on all sample points
if probinfo.analyticflag;
    contgrdsamples = feval(probinfo.contgrd, continput);
else
    contgrdsamples = feval(probinfo.contgrd, continput, probinfo);
end
% -------------------------------------------------------------------------

% evaluate endpoint function at random points -----------------------------
% preallocate the input structure for endpoint functions
ephase(numphase).initialstate = [];
ephase(numphase).finalstate = [];
ephase(numphase).initialtime = [];
ephase(numphase).finaltime = [];
if sum(numintegral,2) ~= 0;
    ephase(numphase).integral = [];
end
% add auxiliary data to endpinput structure
if probinfo.auxflag;
    endpinput.auxdata = probinfo.auxdata;
end

% preallocate derivative output for all sample points
objgrdsamples = zeros(numsamples,probinfo.derivativemap.objnnz1);
if numeventgroup ~= 0;
    eventnnz1 = probinfo.derivativemap.eventnnz1;
    eventgroup(numeventgroup).eventgrd = [];
    for eventgroupcount = 1:numeventgroup;
        eventgroup(eventgroupcount).eventgrd = zeros(numsamples,eventnnz1(eventgroupcount));
    end
end

% evaluate endpoint function derivative at each sample point
for samplecount = 1:numsamples;
    % get endpoint function input for each sample
    for phasecount = 1:numphase;
        ephase(phasecount).initialstate = ephaserand(phasecount).initialstate(samplecount,:);
        ephase(phasecount).finalstate = ephaserand(phasecount).finalstate(samplecount,:);
        ephase(phasecount).initialtime = ephaserand(phasecount).initialtime(samplecount,:);
        ephase(phasecount).finaltime = ephaserand(phasecount).finaltime(samplecount,:);
        if numintegral(phasecount) ~= 0;
            ephase(phasecount).integral = ephaserand(phasecount).integral(samplecount,:);
        end
    end
    endpinput.phase = ephase;
    if numparameters ~= 0;
        endpinput.parameter = pararand(samplecount,:);
    end
    
    % eval endpoint function derivative
    if probinfo.analyticflag;
        endpfirstder = feval(probinfo.endpgrd, endpinput);
    else
        endpfirstder = feval(probinfo.endpgrd, endpinput, probinfo);
    end
    
    % store derivative values of objective for current sample
    objgrdsamples(samplecount,:) = endpfirstder.objectivegrd;
    
    % store derivative values of events for current sample
    if numeventgroup ~= 0;
        for eventgroupcount = 1:numeventgroup;
            eventgroup(eventgroupcount).eventgrd(samplecount,:) = endpfirstder.eventgroup(eventgroupcount).eventgrd;
        end
    end
end
endpgrdsamples.objectivegrd = objgrdsamples;
if numeventgroup ~= 0;
    endpgrdsamples.eventgroup = eventgroup;
end