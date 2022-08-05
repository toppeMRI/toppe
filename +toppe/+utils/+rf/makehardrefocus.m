function [rf, g] = makehardrefocus(flip, dur, nCyclesSpoil, voxWidth, maxSlew, maxGrad)
% function [rf, g] = makehardrefocus(flip, dur, nCyclesSpoil, voxWidth, maxSlew, maxGrad)
%
% Make hard pulse surrounded by spoiler gradients
%
% Inputs
%   flip            [1 1]    degrees
%   dur             [1 1]    sec
%   nCyclesSpoil    [1 1]    number of cycles of gradient crushing across voxWidth
%   voxWidth        [1 1]    voxel size (cm)
%   maxSlew         [1 1]    max slew on one gradient axis
%   maxGrad         [1 1]    max gradient amplitude on one gradient axis
%
% Outputs
%   rf              [nt 1]   complex rf waveform (Gauss)
%   g               [g 1]    gradient waveform (Gauss/cm)

rfhard = toppe.utils.rf.makehardpulse(flip, dur);

gspoil = toppe.utils.makecrusher(nCyclesSpoil, voxWidth, 0, maxSlew/sqrt(2), maxGrad/sqrt(2));
gspoil = [0; 0; gspoil(:); 0; 0];  % to reduce change of overlap between gradients and rf

rf = [0*gspoil(:);   rfhard(:); 0*gspoil(:)];
g  = [  gspoil(:); 0*rfhard(:);   gspoil(:)];
rf = toppe.utils.makeGElength(rf);
g  = toppe.utils.makeGElength(g);
% toppe.writemod('rf', rf, 'gx', g, 'gy', g, 'ofname', mods.refocus);

