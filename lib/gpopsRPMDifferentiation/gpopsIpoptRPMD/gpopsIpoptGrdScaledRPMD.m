function Grdout = gpopsIpoptGrdScaledRPMD(Zin, probinfo)

% unscale variables
Zunscaled = (Zin - probinfo.Zshift)./probinfo.Zscale;

% get gradient nonzeros
grdnz = probinfo.grdscale.*gpopsGrdnzRPMD(Zunscaled, probinfo);

% full gradient
Grdout = zeros(1,probinfo.nlpnumvar);
Grdout(probinfo.grdpat) = grdnz;
%}