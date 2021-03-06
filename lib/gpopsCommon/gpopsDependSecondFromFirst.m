function derivativemap = gpopsDependSecondFromFirst(derivativemap, probinfo)

% gpopsDependSecondFromFirst
% this function creates the second derivative map from the first derivative
% map. The nonzero locations of the second derivative map are estimated
% from the first derivative nonzero locations.

% get OCP sizes
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numintegral = probinfo.numintegral;
numpath = probinfo.numpath;
numeventgroup = probinfo.numeventgroup;
numevent = probinfo.numevent;

% preallocate contmap for second derivatives
contnvc2 = zeros(1,numphase);
contmap2(numphase).contvarmap2 = [];
contmap2(numphase).dynamicsmap2 = [];
dynamicsnnz2 = zeros(1,numphase);
if sum(numpath,2) ~= 0;
    contmap2(numphase).pathmap2 = [];
    pathnnz2 = zeros(1,numphase);
end
if sum(numintegral,2) ~= 0;
    contmap2(numphase).integrandmap2 = [];
    integrandnnz2 = zeros(1,numphase);
end

% find continuous function second derivative map
for phasecount = 1:numphase;
    % phase information
    numstatep = numstate(phasecount);
    numpathp = numpath(phasecount);
    numintegralp = numintegral(phasecount);
    
    % get first derivative map for phase
    contmap1 = derivativemap.contmap1(phasecount);
    
    % get dynamics first sparsity
    dynamicssp1 = double(contmap1.dynamicsmap1 ~= 0);
    
    % contMsp2 is the second derivative sparsity pattern of all functions for
    % phase
    contMsp2 = tril(dynamicssp1'*dynamicssp1);
    
    % get path first sparsity
    if numpathp ~= 0;
        pathsp1 = double(contmap1.pathmap1 ~= 0);
        
        % add path to contMsp2
        contMsp2 = contMsp2 + tril(pathsp1'*pathsp1);
    end
    
    % get integrand first sparsity
    if numintegralp ~= 0;
        intsp1 = double(contmap1.integrandmap1 ~= 0);
        
        % add integrand to contMsp2
        contMsp2 = contMsp2 + tril(intsp1'*intsp1);
    end
    
    % find second derivative variable map
    [row, col] = find(contMsp2);
    contlocmap2 = sortrows([row, col])';
    
    % get proper variable numbers
    contmap2(phasecount).contvarmap2 = contmap1.contvarmap1(contlocmap2);
    
    % get number of partials to be taken of the continuous function to find
    % all nonzero second derivatives
    contnvc2(phasecount) = size(contlocmap2,2);
    
    % reference used to find dynamics, path and integrand second derivative maps
    contsp2ref = sub2ind(size(contMsp2), contlocmap2(1,:), contlocmap2(2,:));
    
    % find second derivative dynamics map
    dynamicsmap2 = zeros(numstatep, contnvc2(phasecount));
    for statecount = 1:numstatep;
        % get second derivative sparsity of dynamics
        dynamicssp2Mat = tril(dynamicssp1(statecount,:)'*dynamicssp1(statecount,:));
        dynamicsmap2(statecount,:) = dynamicssp2Mat(contsp2ref);
    end
    % number of dynamics nonzero second derivatives
    dynamicsnnz2(phasecount) = nnz(dynamicsmap2);
    
    % get dynamics second derivative locations
    dynamicsmap2(dynamicsmap2 ~= 0) = 1:dynamicsnnz2(phasecount);
    
    % save dynamicsmap2 in structure
    contmap2(phasecount).dynamicsmap2 = dynamicsmap2;
    
    % find second derivative path map
    if numpathp ~= 0;
        pathmap2 = zeros(numpathp, contnvc2(phasecount));
        for pathcount = 1:numpathp;
            % get second derivative sparsity of path
            pathsp2Mat = tril(pathsp1(pathcount,:)'*pathsp1(pathcount,:));
            pathmap2(pathcount,:) = pathsp2Mat(contsp2ref);
        end
        % number of path nonzero second derivatives
        pathnnz2(phasecount) = nnz(pathmap2);
        
        % get path second derivative locations
        pathmap2(pathmap2 ~= 0) = 1:pathnnz2(phasecount);
        
        % save pathmap2 in structure
        contmap2(phasecount).pathmap2 = pathmap2;
    end
    
    % find second derivative integrand map
    if numintegralp ~= 0;
        integrandmap2 = zeros(numintegralp, contnvc2(phasecount));
        for intcount = 1:numintegralp;
            % get second derivative sparsity of integrand
            intsp2Mat = tril(intsp1(intcount,:)'*intsp1(intcount,:));
            integrandmap2(intcount,:) = intsp2Mat(contsp2ref);
        end
        % number of integrand nonzero second derivatives
        integrandnnz2(phasecount) = nnz(integrandmap2);
        
        % get integrand second derivative locations
        integrandmap2(integrandmap2 ~= 0) = 1:integrandnnz2(phasecount);
        
        % save integrandmap2 in structure
        contmap2(phasecount).integrandmap2 = integrandmap2;
    end
end

% endpoint second derivative map
% get objective first sparsity
objsp1 = double(derivativemap.endpobjmap1 ~= 0);

% endpMsp2 is the second derivative sparsity pattern of all functions of
% the endpoint function
endpMsp2 = tril(objsp1'*objsp1);

% get event first sparsity
if numeventgroup ~= 0;
    % Get first event map
    endpeventmap1 = derivativemap.endpeventmap1;
    
    eventgroupsp2(numeventgroup).eventsp1 = [];
    % event second derivative function map
    for eventgroupcount = 1:numeventgroup;
        if numevent(eventgroupcount) ~= 0;
            % Get first event map
            eventsp1 = double(endpeventmap1(eventgroupcount).first ~= 0);
            
            % add event to endpMsp2
            endpMsp2 = endpMsp2 + tril(eventsp1'*eventsp1);
            
            % store eventsp1 in structure
            eventgroupsp2(eventgroupcount).eventsp1 = eventsp1;
        end
    end
end

% find endp second derivative variable map
[row, col] = find(endpMsp2);
endplocmap2 = sortrows([row, col])';

% get proper variable numbers
endpvarmap2 = derivativemap.endpvarmap1(endplocmap2);

% get number of partials to be taken of the endpoint function to find
% all nonzero second derivatives
endpnvc2 = size(endplocmap2,2);

% reference used to find event second derivative map
endpsp2ref = sub2ind(size(endpMsp2), endplocmap2(1,:), endplocmap2(2,:));


% find second derivative objective map
% get second derivative sparsity of objective
objsp2Mat = tril(objsp1'*objsp1);
objfunmap2 = objsp2Mat(endpsp2ref);

% number of objective nonzero second derivatives
objnnz2 = nnz(objfunmap2);

% get objective second derivative locations
objfunmap2(objfunmap2 ~= 0) = 1:objnnz2;

% event second derivative map
if numeventgroup ~= 0;    
    % preallocate eventfunmap2
    endpeventmap2(numeventgroup).second = [];
    eventnnz2 = zeros(1,numeventgroup);
    
    % event second derivative function map
    for eventgroupcount = 1:numeventgroup;
        if numevent(eventgroupcount) ~= 0;
            eventmap2 = zeros(numevent(eventgroupcount), endpnvc2);
            for eventcount = 1:numevent(eventgroupcount);
                % get second derivative sparsity of event
                eventsp1 = eventgroupsp2(eventgroupcount).eventsp1(eventcount,:);
                eventsp2Mat = tril(eventsp1'*eventsp1);
                eventmap2(eventcount,:) = eventsp2Mat(endpsp2ref);
            end
            % number of event nonzero second derivatives
            eventnnz2(eventgroupcount) = nnz(eventmap2);
            
            % get event second derivative locations
            eventmap2(eventmap2 ~= 0) = 1:eventnnz2(eventgroupcount);
            
            % save eventmap2 in structure
            endpeventmap2(eventgroupcount).second = eventmap2;            
        end
    end
end

% second derivative map
% second derivative endp map
derivativemap.endpnvc2 = endpnvc2;
derivativemap.endpvarmap2 = endpvarmap2;
derivativemap.endpobjmap2 = objfunmap2;
if numeventgroup ~= 0;
    derivativemap.endpeventmap2 = endpeventmap2;
end

% second derivative obj and event nnz
derivativemap.objnnz2 = objnnz2;
if numeventgroup ~= 0;
    derivativemap.eventnnz2 = eventnnz2;
end

% second derivative cont map
derivativemap.contnvc2 = contnvc2;
derivativemap.contmap2 = contmap2;

% second derivative dynamics, path, and intergrand nnz
derivativemap.dynamicsnnz2 = dynamicsnnz2;
if sum(numpath,2) ~= 0;
    derivativemap.pathnnz2 = pathnnz2;
end
if sum(numintegral,2) ~= 0;
    derivativemap.integrandnnz2 = integrandnnz2;
end