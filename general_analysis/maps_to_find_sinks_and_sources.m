function [ld_rg_counts_orig, ld_yv_counts_orig, yv_rg_counts_orig, ld_rg_counts_filter, ld_yv_counts_filter, yv_rg_counts_filter] = maps_to_find_sinks_and_sources(orig_data, filter_data)
ldmino = min(orig_data(:, 1)) - 0.02;
rgmino = min(orig_data(:, 2)) - 0.02;
yvmino = min(orig_data(:, 3)) - 0.02;

ldmaxo = max(orig_data(:, 1)) + 0.02;
rgmaxo = max(orig_data(:, 2)) + 0.02;
yvmaxo = max(orig_data(:, 3)) + 0.02;

ldminf = min(filter_data(:, 1)) - 0.02;
rgminf = min(filter_data(:, 2)) - 0.02;
yvminf = min(filter_data(:, 3)) - 0.02;

ldmaxf = max(filter_data(:, 1)) + 0.02;
rgmaxf = max(filter_data(:, 2)) + 0.02;
yvmaxf = max(filter_data(:, 3)) + 0.02;

ldmin = min(ldmino, ldminf);
rgmin = min(rgmino, rgminf);
yvmin = min(yvmino, yvminf);

ldmax = max(ldmaxo, ldmaxf);
rgmax = max(rgmaxo, rgmaxf);
yvmax = max(yvmaxo, yvmaxf);

steps = 300;
ld_sps = linspace(ldmin, ldmax, steps);
rg_sps = linspace(rgmin, rgmax, steps);
yv_sps = linspace(yvmin, yvmax, steps);
ld_rg_counts_orig = zeros(steps, steps);
ld_yv_counts_orig = zeros(steps, steps);
yv_rg_counts_orig = zeros(steps, steps);
ld_rg_counts_filter = zeros(steps, steps);
ld_yv_counts_filter = zeros(steps, steps);
yv_rg_counts_filter = zeros(steps, steps);
lds = squeeze(orig_data(:, 1));
rgs = squeeze(orig_data(:, 2));
yvs = squeeze(orig_data(:, 3));
step_size_ld = abs(ld_sps(2) - ld_sps(1));
step_size_rg = abs(rg_sps(2) - rg_sps(1));
step_size_yv = abs(yv_sps(2) - yv_sps(1));

% do it first for the original data
for ldc = 1:length(ld_sps)
    ld = ld_sps(ldc) + step_size_ld/2;
    for rgc = 1:length(rg_sps)
        rg = rg_sps(rgc) + step_size_rg/2;
        d1_idxs = ((lds > ld - step_size_ld/2) & (lds < ld + step_size_ld/2)) & ((rgs > rg - step_size_rg/2) & (rgs < rg + step_size_rg/2));
        
        if ~isempty(find(d1_idxs))
            ld_rg_counts_orig(ldc, rgc) = length(d1_idxs);
        end
    end
end

for yvc = 1:length(yv_sps)
    yv = yv_sps(yvc) + step_size_yv/2;
    for rgc = 1:length(rg_sps)
        rg = rg_sps(rgc) + step_size_rg/2;
        d2_idxs = ((rgs > rg - step_size_rg/2) & (rgs < rg + step_size_rg/2)) & ((yvs > yv - step_size_yv/2) & (yvs < yv + step_size_yv/2));
        
        if ~isempty(find(d2_idxs))
            yv_rg_counts_orig(yvc, rgc) = length(d2_idxs);
        end
    end
end

for ldc = 1:length(ld_sps)
    ld = ld_sps(ldc) + step_size_ld/2;
    for yvc = 1:length(yv_sps)
        yv = yv_sps(yvc) + step_size_yv/2;
        d3_idxs = ((lds > ld - step_size_ld/2) & (lds < ld + step_size_ld/2)) & ((yvs > yv - step_size_yv/2) & (yvs < yv + step_size_yv/2));
        
        if ~isempty(find(d3_idxs))
            ld_yv_counts_orig(ldc, yvc) = length(d3_idxs);
        end
    end
end

% then, do it for the filter data
lds = squeeze(filter_data(:, 1));
rgs = squeeze(filter_data(:, 2));
yvs = squeeze(filter_data(:, 3));

for ldc = 1:length(ld_sps)
    ld = ld_sps(ldc) + step_size_ld/2;
    for rgc = 1:length(rg_sps)
        rg = rg_sps(rgc) + step_size_rg/2;
        d1_idxs = ((lds > ld - step_size_ld/2) & (lds < ld + step_size_ld/2)) & ((rgs > rg - step_size_rg/2) & (rgs < rg + step_size_rg/2));
        
        if ~isempty(find(d1_idxs))
            ld_rg_counts_filter(ldc, rgc) = length(d1_idxs);
        end
    end
end

for yvc = 1:length(yv_sps)
    yv = yv_sps(yvc) + step_size_yv/2;
    for rgc = 1:length(rg_sps)
        rg = rg_sps(rgc) + step_size_rg/2;
        d2_idxs = ((rgs > rg - step_size_rg/2) & (rgs < rg + step_size_rg/2)) & ((yvs > yv - step_size_yv/2) & (yvs < yv + step_size_yv/2));
        
        if ~isempty(find(d2_idxs))
            yv_rg_counts_filter(yvc, rgc) = length(d2_idxs);
        end
    end
end

for ldc = 1:length(ld_sps)
    ld = ld_sps(ldc) + step_size_ld/2;
    for yvc = 1:length(yv_sps)
        yv = yv_sps(yvc) + step_size_yv/2;
        d3_idxs = ((lds > ld - step_size_ld/2) & (lds < ld + step_size_ld/2)) & ((yvs > yv - step_size_yv/2) & (yvs < yv + step_size_yv/2));
        
        if ~isempty(find(d3_idxs))
            ld_yv_counts_filter(ldc, yvc) = length(d3_idxs);
        end
    end
end

return;
end
