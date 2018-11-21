function w = mdcf_voronoi(kspace, ratio)
%function w = mdcf_voronoi(kspace)
%Inputs:
% - kspace (nk, nd) coordinates that describes the k-space sampling points,
%    the coordinates needs to be normalized to between [-0.5, 0.5]
%    Notice, only input the kspace part that its corresponding acquired
%    signal will be used into reconstruction
% - ratio (1,), dflt: 0.1; threshold the largest $ratio of w before outputing.
%Outputs:
% - w (nk,) dcf output
%

import toppe.utils.spiral.*

if nargin == 0, w = test(); return; end
if nargin == 1, ratio = 0.1; end
if any(abs(kspace)>0.5)
  dims = round(max(kspace,[],1)) - round(min(kspace,[],1)) + 1;
  kspace = bsxfun(@rdivide, kspace, dims);
  warning('kspace normalized to interval [-0.5, 0.5]');
end

%% calc w
[kspaceU, ~, ik] = unique(kspace, 'rows', 'stable');

[v, c] = voronoin(kspaceU);
nk = size(kspaceU, 1);

wi = zeros(nk, 1);

getV = @(inds)v(inds,:);
x = cellfun(getV, c, 'UniformOutput', false);

parforWarning();
parfor ii = 1:nk % parfor ii = 1:nk
  try [~, wi(ii)] = convhulln(x{ii});
  catch % x{ii} gen'ed w/ voronoi thus must be convex, here deals inf case
    wi(ii) = inf;
  end
end

wo = wi(ik);
[bincount, ind] = histc(ik, unique(ik));
wo = wo./bincount(ind);

%% clean up
a = 0.98;
while true
  r = (1 - (1-ratio)*a);
  thld = selectNth(wo, round(r * numel(wo))); % 0.1 chosen arbitrarily
  if ~isfinite(thld), a = a*a;
  else, break, end
end
wo(wo>thld) = thld;

wo(isnan(wo)) = 0;
w = wo ./ max(wo);
end

function w = test()

load kSpiral2D k % load kSpiral2D imSize
k = [real(k(:)), imag(k(:))]; % already normalized
% k_cycleFov = bsxfun(@times, k, imSize);
w = mdcf_voronoi(k, 0.1);

end
