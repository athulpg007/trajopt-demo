function contsecondder = gpopsContSecondFD(input, probinfo)

% gpopsContSecondFD
% this function estimates the second derivatives of the OCP continuous
% function using forward difference
% only the derivatives indicated in the derivativemap are found
% uses value based step sizes

% get OCP info
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numintegral = probinfo.numintegral;
numpath = probinfo.numpath;
pathsum = sum(numpath,2);
integrandsum = sum(numintegral,2);
numparameters = probinfo.numparameter;

% get continuous function derivative map
contmap2 = probinfo.derivativemap.contmap2;

% get nvc for second derivatives
contnvc = probinfo.derivativemap.contnvc2;
maxnvc = max(contnvc);

% preallocate derivative output
% only nonzero derivatives of the optimal control problem are stored
dynamicsnnz = probinfo.derivativemap.dynamicsnnz2;
contsecondder(numphase).dynamicshes = [];
if pathsum ~= 0;
    pathnnz = probinfo.derivativemap.pathnnz2;
    contsecondder(numphase).pathhes = [];
else
    pathnnz = zeros(1,numphase);
end
if integrandsum ~= 0;
    integrandnnz = probinfo.derivativemap.integrandnnz2;
    contsecondder(numphase).integrandhes = [];
else
    integrandnnz = zeros(1,numphase);
end
for phasecount = 1:numphase;
    NN = size(input.phase(phasecount).time,1);
    if dynamicsnnz(phasecount) ~= 0;
        contsecondder(phasecount).dynamicshes = zeros(NN,dynamicsnnz(phasecount));
    end
    if pathnnz(phasecount) ~= 0;
        contsecondder(phasecount).pathhes = zeros(NN,pathnnz(phasecount));
    end
    if integrandnnz(phasecount) ~= 0;
        contsecondder(phasecount).integrandhes = zeros(NN,integrandnnz(phasecount));
    end
end

% get base stepsize
ustep = probinfo.stepsize;

% preallocate hstep
hstep(numphase).h1 = [];
hstep(numphase).h2 = [];

% get non-perturbed function solution
output = feval(probinfo.contfunction, input);

% find derivatives of all phases simultaneously
for nvccount = 1:maxnvc;
    % initiate pertinput as the unperturbed input
    pertinput1 = input;
    pertinput2 = input;
    pertinput3 = input;
    for phasecount = 1:numphase;
        numstatep = numstate(phasecount);
        numcontrolp = numcontrol(phasecount);
        % test if the perturbed variable count is less then the number of
        % perturbed variables in the phase
        if nvccount <= contnvc(1,phasecount);
            varnum1 = contmap2(phasecount).contvarmap2(1,nvccount);
            varnum2 = contmap2(phasecount).contvarmap2(2,nvccount);
            % perturb first variable
            if varnum1 <= numstatep;
                % perturb state
                refmark = varnum1;
                h = ustep.*(abs(input.phase(phasecount).state(:,refmark))+1);
                hstep(phasecount).h1 = h;
                pertinput1.phase(phasecount).state(:,refmark) = pertinput1.phase(phasecount).state(:,refmark) + h;
                pertinput2.phase(phasecount).state(:,refmark) = pertinput2.phase(phasecount).state(:,refmark) + h;
            elseif varnum1 <= numstatep+numcontrolp;
                % perturb control
                refmark = varnum1-numstatep;
                h = ustep.*(abs(input.phase(phasecount).control(:,refmark))+1);
                hstep(phasecount).h1 = h;
                pertinput1.phase(phasecount).control(:,refmark) = pertinput1.phase(phasecount).control(:,refmark) + h;
                pertinput2.phase(phasecount).control(:,refmark) = pertinput2.phase(phasecount).control(:,refmark) + h;
            elseif varnum1 == numstatep+numcontrolp+1;
                % perturb time
                h = ustep.*(abs(input.phase(phasecount).time)+1);
                hstep(phasecount).h1 = h;
                pertinput1.phase(phasecount).time = pertinput1.phase(phasecount).time + h;
                pertinput2.phase(phasecount).time = pertinput2.phase(phasecount).time + h;
            elseif varnum1 <= numstatep+numcontrolp+1+numparameters;
                % perturb parameter
                refmark = varnum1-numstatep-numcontrolp-1;
                h = ustep.*(abs(input.phase(phasecount).parameter(:,refmark))+1);
                hstep(phasecount).h1 = h;
                pertinput1.phase(phasecount).parameter(:,refmark) = pertinput1.phase(phasecount).parameter(:,refmark) + h;
                pertinput2.phase(phasecount).parameter(:,refmark) = pertinput2.phase(phasecount).parameter(:,refmark) + h;
            end
            % perturb second variable
            if varnum2 <= numstatep;
                % perturb state
                refmark = varnum2;
                h = ustep.*(abs(input.phase(phasecount).state(:,refmark))+1);
                hstep(phasecount).h2 = h;
                pertinput1.phase(phasecount).state(:,refmark) = pertinput1.phase(phasecount).state(:,refmark) + h;
                pertinput3.phase(phasecount).state(:,refmark) = pertinput3.phase(phasecount).state(:,refmark) + h;
            elseif varnum2 <= numstatep+numcontrolp;
                % perturb control
                refmark = varnum2-numstatep;
                h = ustep.*(abs(input.phase(phasecount).control(:,refmark))+1);
                hstep(phasecount).h2 = h;
                pertinput1.phase(phasecount).control(:,refmark) = pertinput1.phase(phasecount).control(:,refmark) + h;
                pertinput3.phase(phasecount).control(:,refmark) = pertinput3.phase(phasecount).control(:,refmark) + h;
            elseif varnum2 == numstatep+numcontrolp+1;
                % perturb time
                h = ustep.*(abs(input.phase(phasecount).time)+1);
                hstep(phasecount).h2 = h;
                pertinput1.phase(phasecount).time = pertinput1.phase(phasecount).time + h;
                pertinput3.phase(phasecount).time = pertinput3.phase(phasecount).time + h;
            elseif varnum2 <= numstatep+numcontrolp+1+numparameters;
                % perturb parameter
                refmark = varnum2-numstatep-numcontrolp-1;
                h = ustep.*(abs(input.phase(phasecount).parameter(:,refmark))+1);
                hstep(phasecount).h2 = h;
                pertinput1.phase(phasecount).parameter(:,refmark) = pertinput1.phase(phasecount).parameter(:,refmark) + h;
                pertinput3.phase(phasecount).parameter(:,refmark) = pertinput3.phase(phasecount).parameter(:,refmark) + h;
            end
        end
    end
    
    % evaluate function on perturbed input
    pertoutput1 = feval(probinfo.contfunction, pertinput1);
    pertoutput2 = feval(probinfo.contfunction, pertinput2);
    pertoutput3 = feval(probinfo.contfunction, pertinput3);

    
    % calculate the derivative value in each phase
    for phasecount = 1:numphase;
        numstatep = numstate(phasecount);
        numpathp = numpath(phasecount);
        numintegralp = numintegral(phasecount);
        % test if the perturbed variable count is less then the number of
        % perturbed variables in the phase
        if nvccount <= contnvc(1,phasecount);
            % calculate nonzero derivatives of dynamic constraints
            for dynamicscount = 1:numstatep;
                refmark = contmap2(phasecount).dynamicsmap2(dynamicscount,nvccount);
                if refmark ~= 0;
                    contsecondder(phasecount).dynamicshes(:,refmark) = (pertoutput1(phasecount).dynamics(:,dynamicscount) - pertoutput2(phasecount).dynamics(:,dynamicscount) - pertoutput3(phasecount).dynamics(:,dynamicscount) + output(phasecount).dynamics(:,dynamicscount))./(hstep(phasecount).h1.*hstep(phasecount).h2);
                end
            end
            % calculate nonzero derivatives of path constraints
            for pathcount = 1:numpathp;
                refmark = contmap2(phasecount).pathmap2(pathcount,nvccount);
                if refmark ~= 0;
                    contsecondder(phasecount).pathhes(:,refmark) = (pertoutput1(phasecount).path(:,pathcount) - pertoutput2(phasecount).path(:,pathcount) - pertoutput3(phasecount).path(:,pathcount) + output(phasecount).path(:,pathcount))./(hstep(phasecount).h1.*hstep(phasecount).h2);
                end
            end
            % calculate nonzero derivatives of intergrand constraints
            for intergralcount = 1:numintegralp;
                refmark = contmap2(phasecount).integrandmap2(intergralcount,nvccount);
                if refmark ~= 0;
                    contsecondder(phasecount).integrandhes(:,refmark) = (pertoutput1(phasecount).integrand(:,intergralcount) - pertoutput2(phasecount).integrand(:,intergralcount) - pertoutput3(phasecount).integrand(:,intergralcount) + output(phasecount).integrand(:,intergralcount))./(hstep(phasecount).h1.*hstep(phasecount).h2);
                end
            end
        end
    end
end