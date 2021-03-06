function [hespat, probinfo] = gpopsHesPatRPMI(probinfo)

% gpopsHesPatRPMI
% this function finds the sparsity pattern of the lower triangle of the 
% Hessian, the output hespat contains row and column indices of the Hessian
% nonzeros, the values associated with duplicate indices will be summed
% together to define the value for that element

% NLP variable order
% [states*(nodes+1); controls*nodes; t0; tf, Q] for each phase
% [stack all phases(1;...;numphase); parameters]

% get OCP info
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numparameter = probinfo.numparameter;

% get parametermap
if numparameter ~= 0;
    parametermap = probinfo.nlpparametermap;
end

% get number of nodes
numnodes = probinfo.numnodes;

% get hesnvc
hesnvc = probinfo.derivativemap.hesnvc;

% preallocate hespat
hespat = zeros(1,2);

% Hessian nonzero locations from continuous function
hesmarkere = 0;
for phasecount = 1:numphase;
    % OCP info for phase
    numstatep = numstate(phasecount);
    numcontrolp = numcontrol(phasecount);
    numnodesp = numnodes(phasecount);
    
    % get NLP contmap for each phase, and derivative map for each phase
    nlpcontmap = probinfo.nlpcontmap(phasecount);
    contmap2 = probinfo.derivativemap.contmap2(phasecount);
    
    % find nonzero locations for each OCP variable
    for contnvc2count = 1:hesnvc(phasecount);
        varnum1 = contmap2.conthesmap(1,contnvc2count);
        varnum2 = contmap2.conthesmap(2,contnvc2count);
        if varnum1 <= numstatep;
            % varnum1 is state
            % varnum2 is state (must be <= varnum1)
            stateref1 = varnum1;
            stateref2 = varnum2;
            % Hessian nonzero locations
            rows = nlpcontmap.statemap(1,stateref1):(nlpcontmap.statemap(2,stateref1)-1);
            cols = nlpcontmap.statemap(1,stateref2):(nlpcontmap.statemap(2,stateref2)-1);
            % get assignment index
            hesmarkers = hesmarkere + 1;
            hesmarkere = hesmarkere + numnodesp;
            hesindex = hesmarkers:hesmarkere;
            % assign values
            hespat(hesindex,1) = rows;
            hespat(hesindex,2) = cols;
        elseif varnum1 <= numstatep+numcontrolp;
            % varnum1 is control
            controlref1 = varnum1-numstatep;
            rows = nlpcontmap.controlmap(1,controlref1):nlpcontmap.controlmap(2,controlref1);
            if varnum2 <= numstatep;
                % varnum2 is a state
                stateref2 = varnum2;
                cols = nlpcontmap.statemap(1,stateref2):(nlpcontmap.statemap(2,stateref2)-1);
            else
                % varnum2 is a control
                controlref2 = varnum2-numstatep;
                cols = nlpcontmap.controlmap(1,controlref2):nlpcontmap.controlmap(2,controlref2);
            end
            % get assignment index
            hesmarkers = hesmarkere + 1;
            hesmarkere = hesmarkere + numnodesp;
            hesindex = hesmarkers:hesmarkere;
            % assign values
            hespat(hesindex,1) = rows;
            hespat(hesindex,2) = cols;
        elseif varnum1 == numstatep+numcontrolp+1;
            % varnum1 is time
            rowt0 = nlpcontmap.timemap(1);
            rowtf = nlpcontmap.timemap(2);
            if varnum2 <= numstatep+numcontrolp;
                if varnum2 <= numstatep;
                    % varnum2 is state
                    stateref2 = varnum2;
                    cols = nlpcontmap.statemap(1,stateref2):(nlpcontmap.statemap(2,stateref2)-1);
                else
                    % varnum2 is control
                    controlref2 = varnum2-numstatep;
                    cols = nlpcontmap.controlmap(1,controlref2):nlpcontmap.controlmap(2,controlref2);
                end
                % get assignment index
                hesmarkers = hesmarkere + 1;
                hesmarkere = hesmarkere + numnodesp;
                hesindex = hesmarkers:hesmarkere;
                % assign values
                hespat(hesindex,1) = rowt0;
                hespat(hesindex,2) = cols;
                % get assignment index
                hesmarkers = hesmarkere + 1;
                hesmarkere = hesmarkere + numnodesp;
                hesindex = hesmarkers:hesmarkere;
                % assign values
                hespat(hesindex,1) = rowtf;
                hespat(hesindex,2) = cols;
            else
                % t0/t0
                hesmarkere = hesmarkere + 1;
                hespat(hesmarkere,1) = nlpcontmap.timemap(1);
                hespat(hesmarkere,2) = nlpcontmap.timemap(1);
                % tf/t0
                hesmarkere = hesmarkere + 1;
                hespat(hesmarkere,1) = nlpcontmap.timemap(2);
                hespat(hesmarkere,2) = nlpcontmap.timemap(1);
                % tf/tf
                hesmarkere = hesmarkere + 1;
                hespat(hesmarkere,1) = nlpcontmap.timemap(2);
                hespat(hesmarkere,2) = nlpcontmap.timemap(2);
            end
        elseif varnum1 <= numstatep+numcontrolp+1+numparameter;
            % varnum1 is parameter
            parameterref1 = varnum1-numstatep-numcontrolp-1;
            rowpara = parametermap(parameterref1);
            if varnum2 <= numstatep+numcontrolp;
                if varnum2 <= numstatep;
                    % varnum2 is state
                    stateref2 = varnum2;
                    cols = nlpcontmap.statemap(1,stateref2):(nlpcontmap.statemap(2,stateref2)-1);
                else
                    % varnum2 is control
                    controlref2 = varnum2-numstatep;
                    cols = nlpcontmap.controlmap(1,controlref2):nlpcontmap.controlmap(2,controlref2);
                end
                % get assignment index
                hesmarkers = hesmarkere + 1;
                hesmarkere = hesmarkere + numnodesp;
                hesindex = hesmarkers:hesmarkere;
                % assign values
                hespat(hesindex,1) = rowpara;
                hespat(hesindex,2) = cols;
            elseif varnum2 == numstatep+numcontrolp+1;
                % varnum2 is time
                % parameter/t0
                hesmarkere = hesmarkere + 1;
                hespat(hesmarkere,1) = rowpara;
                hespat(hesmarkere,2) = nlpcontmap.timemap(1);
                % parameter/tf
                hesmarkere = hesmarkere + 1;
                hespat(hesmarkere,1) = rowpara;
                hespat(hesmarkere,2) = nlpcontmap.timemap(2);
            else
                % varnum2 is parameter
                parameterref2 = varnum2-numstatep-numcontrolp-1;
                colpara = parametermap(parameterref2);
                % parameterI/parameterJ
                hesmarkere = hesmarkere + 1;
                hespat(hesmarkere,1) = rowpara;
                hespat(hesmarkere,2) = colpara;
            end
        end
    end
end

% get endpnvc2, nlpendpvarmap and endpoint second derivative map
endpnvc2 = probinfo.derivativemap.endpnvc2;
nlpendpvarmap = probinfo.nlpendpvarmap;
endpvarmap2 = probinfo.derivativemap.endpvarmap2;

% Hessian nonzero locations from endpoint function
for endpnvccount = 1:endpnvc2;
    varnum1 = endpvarmap2(1,endpnvccount);
    varnum2 = endpvarmap2(2,endpnvccount);
    % assign location
    nlpvar1 = nlpendpvarmap(varnum1);
    nlpvar2 = nlpendpvarmap(varnum2);
    hesmarkere = hesmarkere + 1;
    if nlpvar1 >= nlpvar2;
        hespat(hesmarkere,1) = nlpvar1;
        hespat(hesmarkere,2) = nlpvar2;
    else
        hespat(hesmarkere,1) = nlpvar2;
        hespat(hesmarkere,2) = nlpvar1;
    end
end

probinfo.hesnnz = size(hespat,1);