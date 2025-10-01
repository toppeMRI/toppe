function view_kdata(kx,ky, kdata, n_cols, title_str)
% Display (multi-shot) stack-of-spiral k data
% Inputs:
%   kx/ky   : [nread,nleaf]
%   kdata   : [nread,nleaf, n_kz]
%   n_cols  : number of columns in the figure, dflt = 4;
%   title_str: (optional) title of figure
%
% Example:
%   view_kdata(kx,ky, kdata, 4, 'SOSP k-Space')

if nargin < 4, n_cols = 4; end
if nargin < 5, title_str = 'Stack of Spiral k Space'; end


n_slices_to_show = size(kdata,3);
n_rows = ceil(n_slices_to_show / n_cols);

% ---- Step 1: Precompute global min and max for color scale ----
mags = abs(kdata(:));
vmin = min(mags);
vmax = max(mags);

%figure('Position',[100 100 300*n_cols 300*n_rows])
figure
for iz = 1:n_slices_to_show
    subplot(n_rows, n_cols, iz);
    
    for ileaf=1:size(kdata,2)
        kd = squeeze(kdata(:,ileaf,iz));

        mag = abs(kd(:));

        % Scatter plot
        scatter(kx(:,ileaf), ky(:,ileaf), 20, mag, 'filled')
    end
    colormap('gray')
    caxis([vmin vmax])
    axis equal
    % xlim([-pi pi])
    % ylim([-pi pi])
    xlabel('kx')
    ylabel('ky')
    %title(sprintf('Slice %d', z_idx(iz)))
end

% Shared colorbar
h = colorbar('Position', [0.92 0.15 0.015 0.7]);
ylabel(h,'Log Magnitude')

sgtitle(title_str, 'FontSize', 16)
end
