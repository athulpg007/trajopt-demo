
% optimalprimePathSetup
% this script adds the appropriate paths for use of optimalprime to the
% MATLAB path directory

% get current directory
currdir = pwd;

% add opRPMintegration/opRPMIsnopt/
disp(['Adding Directory ',pwd,'/lib/gpopsRPMIntegration/gpopsSnoptRPMI/ to Path']);
addpath([pwd,'/lib/gpopsRPMIntegration/gpopsSnoptRPMI/'],'-begin');

% add opRPMintegration/opRPMIipopt/
disp(['Adding Directory ',pwd,'/lib/gpopsRPMIntegration/gpopsIpoptRPMI/ to Path']);
addpath([pwd,'/lib/gpopsRPMIntegration/gpopsIpoptRPMI/'],'-begin');

% add opRPMintegration/
disp(['Adding Directory ',pwd,'/lib/gpopsRPMIntegration/ to Path']);
addpath([pwd,'/lib/gpopsRPMIntegration/'],'-begin');

% add opRPMdifferentiation/opRPMDsnopt/
disp(['Adding Directory ',pwd,'/lib/gpopsRPMDifferentiation/gpopsSnoptRPMD/ to Path']);
addpath([pwd,'/lib/gpopsRPMDifferentiation/gpopsSnoptRPMD/'],'-begin');

% add opRPMdifferentiation/opRPMDipopt/
disp(['Adding Directory ',pwd,'/lib/gpopsRPMDifferentiation/gpopsIpoptRPMD/ to Path']);
addpath([pwd,'/lib/gpopsRPMDifferentiation/gpopsIpoptRPMD/'],'-begin');

% add opRPMdifferentiation/
disp(['Adding Directory ',pwd,'/lib/gpopsRPMDifferentiation/ to Path']);
addpath([pwd,'/lib/gpopsRPMDifferentiation/'],'-begin');
% add opOCPfinitediff/
disp(['Adding Path ',pwd,'/lib/gpopsFiniteDifference/']);
addpath([pwd,'/lib/gpopsFiniteDifference/'],'-begin');

% add opMeshHP1/
disp(['Adding Directory ',pwd,'/lib/gpopsMeshHP1/ to Path']);
addpath([pwd,'/lib/gpopsMeshHP1/'],'-begin');

% add opMeshHP/
disp(['Adding Directory ',pwd,'/lib/gpopsMeshHP/ to Path']);
addpath([pwd,'/lib/gpopsMeshHP/'],'-begin');

% add opCommon/
disp(['Adding Directory ',pwd,'/lib/gpopsCommon/ to Path']);
addpath([pwd,'/lib/gpopsCommon/'],'-begin');

% add NLP solver directory
if isdir('nlp/snopt'),
  disp(['Adding Directory ',pwd,'/nlp/snopt/ to Path']);
  addpath([pwd,'/nlp/snopt/'],'-begin');
end;
if isdir('nlp/ipopt'),
  disp(['Adding Directory ',pwd,'/nlp/ipopt/ to Path']);
  addpath([pwd,'/nlp/ipopt/'],'-begin');
end;

savepath;
