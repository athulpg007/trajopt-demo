function ZG = gpopsGuessRPMD(setup, probinfo)

% gpopsGuessRPMD
% this function generates the guess for the nonlinear program from the
% guess given for the optimal control problem found in the setup

% get OCP info
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numintegral = probinfo.numintegral;
numparameter = probinfo.numparameter;

% preallocate ZG
ZG = zeros(probinfo.nlpnumvar,1);

for phasecount = 1:numphase;
    % get guess for phase
    guessp = setup.guess.phase(phasecount);
    
    % get contmap for phase
    contmap = probinfo.nlpcontmap(phasecount);
    
    % OCP info for phase
    numstatep = numstate(phasecount);
    numcontrolp = numcontrol(phasecount);
    numintegralp = numintegral(phasecount);
    guesslength = size(guessp.time,1);
    
    % time domain for guess
    timeguess = guessp.time;
    t0R = timeguess(1,1);
    tfR = timeguess(guesslength,1);
    s = probinfo.collocation(phasecount).s(:,1);
    sp1 = [s; 1];
    timeinterp = (sp1 + 1).*(tfR - t0R)./2 + t0R;
    
    % get NLP guess for state variables and defect constraints in each phase
    for statecount = 1:numstatep;
        % get index from map
        index = contmap.statemap(1,statecount):contmap.statemap(2,statecount);
        
        % interp state guess
        ZG(index) = interp1(timeguess, guessp.state(:,statecount), timeinterp,'pchip');
    end
    
    % get NLP guess for control variables in each phase
    if numcontrolp ~= 0;
        for controlcount = 1:numcontrolp;
            % get index from map
            index = contmap.controlmap(1,controlcount):contmap.controlmap(2,controlcount);
            
            % interp state guess
            ZG(index) = interp1(timeguess, guessp.control(:,controlcount), timeinterp(1:(end-1)),'pchip');
        end
    end
    
    % get NLP guess for time
    % get index for initial time from map
    index = contmap.timemap(1);
    
    % guess for initial time in phase
    ZG(index) = guessp.time(1);
    
    % get index for final time
    index = contmap.timemap(2);
    
    % guess for final time in phase
    ZG(index) = guessp.time(end);
    
    % get NLP guess for integral variables in each phase
    if numintegralp ~= 0;
        for integralcount = 1:numintegralp;
            % get index from map
            index = contmap.integralvarmap(1,integralcount);
            
            % interp integral guess
            ZG(index) = guessp.integral(integralcount);
        end
    end
end

% get NLP guess for parameters
if numparameter ~= 0;
    % getnlp parameter variable map
    parametermap = probinfo.nlpparametermap;
    parameterguess = setup.guess.parameter;
    for parametercount = 1:numparameter;
        % get index
        index = parametermap(parametercount);
        
        % get guess for parameter
        ZG(index) = parameterguess(parametercount);
    end
end