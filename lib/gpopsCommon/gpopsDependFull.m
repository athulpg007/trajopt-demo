function derivativemap = gpopsDependFull(probinfo, derivativelevel)

% gpopsDependFull
% This function gets the optimal control problem dependencies for either
% the first or second derivative levels assuming all functions of the
% optimal control problem have a derivative with respect to every variable
% in the problem

% note
% ___nvc = number of variables that have derivatives with respect to them

% get OCP sizes
numphase = probinfo.numphase;
numstate = probinfo.numstate;
numcontrol = probinfo.numcontrol;
numintegral = probinfo.numintegral;
numpath = probinfo.numpath;
numparameters = probinfo.numparameter;
numeventgroup = probinfo.numeventgroup;

% total number of endpoint variables
numOCPendpvar = 2*sum(numstate,2)+sum(numintegral,2)+2*numphase+numparameters;

% endpoint first derivative variable map
endpvarmap1 = 1:numOCPendpvar;
endpobjmap1 = endpvarmap1;

% objective first derivative map
objvarmap1 = endpvarmap1;
objfunmap1 = endpvarmap1;
objnnz1 = numOCPendpvar;

% event variable first derivative map
if numeventgroup ~= 0;
    % event first derivative variable map
    eventvarmap1 = endpvarmap1;
    
    % preallocate eventfunmap1
    endpeventmap1(numeventgroup).first = [];
    eventnnz1 = zeros(1,numeventgroup);
    
    % event first derivative function map
    for eventgroupcount = 1:numeventgroup;
        % get number of event constraints
        numevent = probinfo.numevent(eventgroupcount);
        endpeventmap1(eventgroupcount).first = reshape(1:numOCPendpvar*numevent,numevent,numOCPendpvar);
        eventnnz1(1,eventgroupcount) = numOCPendpvar*numevent;
    end
    eventfunmap1 = endpeventmap1;
end

% endpoint second derivative map
if derivativelevel == 2;
    % number of elements in lower triangle
    numOCPendptri = (numOCPendpvar.^2 + numOCPendpvar)/2;
    
    % preallocate endpoint second derivative var map
    endpvarmap2 = zeros(2,numOCPendptri);
    
    % endpoint second derivative variable map
    locend = 0;
    for endvarcount = 1:numOCPendpvar;
        endprow = endvarcount*ones(1,endvarcount);
        endpcol = 1:endvarcount;
        locstart = locend + 1;
        locend   = locstart + endvarcount - 1;
        endpvarmap2(:,locstart:locend) = [endprow; endpcol];
    end
    
    % objective second derivative map
    objfunmap2 = 1:numOCPendptri;
    objnnz2 = numOCPendptri;
    
    % event second derivative map
    if numeventgroup ~= 0;
        % preallocate eventfunmap2
        endpeventmap2(numeventgroup).second = [];
        eventnnz2 = zeros(1,numeventgroup);
        
        % event second derivative function map
        for eventgroupcount = 1:numeventgroup;
            % get number of event constraints
            numevent = probinfo.numevent(eventgroupcount);
            endpeventmap2(eventgroupcount).second = reshape(1:numOCPendptri*numevent,numevent,numOCPendptri);
            eventnnz2(1,eventgroupcount) = numOCPendptri*numevent;
        end
    end
end

% preallocate contmap for first derivatives
contmap1(numphase).contvarmap1 = [];
contmap1(numphase).dynamicsmap1 = [];
dynamicsnnz1 = zeros(1,numphase);
contnvc1 = zeros(1,numphase);
if sum(numpath,2) ~= 0;
    contmap1(numphase).pathmap1 = [];
    pathnnz1 = zeros(1,numphase);
end
if sum(numintegral,2) ~= 0;
    contmap1(numphase).integrandmap1 = [];
    integrandnnz1 = zeros(1,numphase);
end

% preallocate contmap for second derivatives
if derivativelevel == 2;
    contmap2(numphase).contvarmap2 = [];
    contmap2(numphase).dynamicsmap2 = [];
    dynamicsnnz2 = zeros(1,numphase);
    contnvc2 = zeros(1,numphase);
    if sum(numpath,2) ~= 0;
        contmap2(numphase).pathmap2 = [];
        pathnnz2 = zeros(1,numphase);
    end
    if sum(numintegral,2) ~= 0;
        contmap2(numphase).integrandmap2 = [];
        integrandnnz2 = zeros(1,numphase);
    end
end

% continuous function derivative map
for phasecount = 1:numphase;
    % get OCP sizes in each phase
    numstatep = numstate(phasecount);
    numcontrolp = numcontrol(phasecount);
    numintegralp = numintegral(phasecount);
    numpathp = numpath(phasecount);
    
    % total number of continuous variables in phase
    numOCPcontvar = numstatep+numcontrolp+1+numparameters;
    
    % get first derivative continuous variable map in phase
    contmap1(phasecount).contvarmap1 = 1:numOCPcontvar;
    
    % get number of variables with derivatives with respect to them
    contnvc1(1,phasecount) = numOCPcontvar;
    
    % get dynamics first derivative map in phase
    contmap1(phasecount).dynamicsmap1 = reshape(1:numOCPcontvar*numstatep,numstatep,numOCPcontvar);
    dynamicsnnz1(1,phasecount) = numOCPcontvar*numstatep;
    
    % find if path constraints are in phase
    if numpathp ~= 0;
        % get path first derivative map in phase
        contmap1(phasecount).pathmap1 = reshape(1:numOCPcontvar*numpathp,numpathp,numOCPcontvar);
        pathnnz1(1,phasecount) = numOCPcontvar*numpathp;
    end
    
    % find if intergal constraints are in phase
    if numintegralp ~= 0;
        % get integrand first derivative map in phase
        contmap1(phasecount).integrandmap1 = reshape(1:numOCPcontvar*numintegralp,numintegralp,numOCPcontvar);
        integrandnnz1(1,phasecount) = numOCPcontvar*numintegralp;
    end
    
    % continuous second derivative map
    if derivativelevel == 2;
        % number of elements in upper triangle
        numOCPconttri = (numOCPcontvar.^2 + numOCPcontvar)/2;
        
        % preallocate endpoint second derivative var map
        contvarmap2 = zeros(2,numOCPconttri);
        
        % endpoint second derivative variable map
        locend = 0;
        for contvarcount = 1:numOCPcontvar;
            controw = contvarcount*ones(1,contvarcount);
            contcol = 1:contvarcount;
            locstart = locend + 1;
            locend   = locstart + contvarcount - 1;
            contvarmap2(:,locstart:locend) = [controw; contcol];
        end
        contmap2(phasecount).contvarmap2 = contvarmap2;
        
        contnvc2(1,phasecount) = size(contvarmap2,2);
        
        % get dynamics second derivative map in phase
        contmap2(phasecount).dynamicsmap2 = reshape(1:numOCPconttri*numstatep,numstatep,numOCPconttri);
        dynamicsnnz2(1,phasecount) = numOCPconttri*numstatep;
        
        % find if path constraints are in phase
        if numpathp ~= 0;
            % get path second derivative map in phase
            contmap2(phasecount).pathmap2 = reshape(1:numOCPconttri*numpathp,numpathp,numOCPconttri);
            pathnnz2(1,phasecount) = numOCPconttri*numpathp;
        end
        
        % find if intergal constraints are in phase
        if numintegralp ~= 0;
            % get integrand second derivative map in phase
            contmap2(phasecount).integrandmap2 = reshape(1:numOCPconttri*numintegralp,numintegralp,numOCPconttri);
            integrandnnz2(1,phasecount) = numOCPconttri*numintegralp;
        end
    end
end

% endp first derivative map
derivativemap.endpnvc1 = size(endpvarmap1,2);
derivativemap.endpvarmap1 = endpvarmap1;
derivativemap.endpobjmap1 = endpobjmap1;
if numeventgroup ~= 0;
    derivativemap.endpeventmap1 = endpeventmap1;
end

% first derivative obj and event nnz
derivativemap.objnnz1 = objnnz1;
if numeventgroup ~= 0;
    derivativemap.eventnnz1 = eventnnz1;
end

% first derivative obj map
derivativemap.objnvc1 = size(objvarmap1,2);
derivativemap.objvarmap1 = objvarmap1;
derivativemap.objfunmap1 = objfunmap1;

% first derivative event map
if numeventgroup ~= 0;
    derivativemap.eventnvc1 = size(eventvarmap1,2);
    derivativemap.eventvarmap1 = eventvarmap1;
    derivativemap.eventfunmap1 = eventfunmap1;
end

% first derivative cont map
derivativemap.contnvc1 = contnvc1;
derivativemap.contmap1 = contmap1;

% first derivative dynamics, path, and intergrand nnz
derivativemap.dynamicsnnz1 = dynamicsnnz1;
if sum(numpath,2) ~= 0;
    derivativemap.pathnnz1 = pathnnz1;
end
if sum(numintegral,2) ~= 0;
    derivativemap.integrandnnz1 = integrandnnz1;
end

% second derivative map
if derivativelevel == 2;
    % second derivative endp map
    derivativemap.endpnvc2 = size(endpvarmap2,2);
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
end