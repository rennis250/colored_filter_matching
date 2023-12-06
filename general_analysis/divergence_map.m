function [ldgr_d1, ldgr_d3, rggr_d1, rggr_d2, yvgr_d2, yvgr_d3, ld_interp_d1, ld_interp_d3, rg_interp_d1, rg_interp_d2, yv_interp_d2, yv_interp_d3] = divergence_map(orig_data, filter_data)
vecs = filter_data - orig_data;

ldmin = min(orig_data(:, 1)) - 0.02;
rgmin = min(orig_data(:, 2)) - 0.02;
yvmin = min(orig_data(:, 3)) - 0.02;

ldmax = max(orig_data(:, 1)) + 0.02;
rgmax = max(orig_data(:, 2)) + 0.02;
yvmax = max(orig_data(:, 3)) + 0.02;

steps = 300;
ld_sps = linspace(ldmin, ldmax, steps);
rg_sps = linspace(rgmin, rgmax, steps);
yv_sps = linspace(yvmin, yvmax, steps);
[rggr_d1, ldgr_d1] = meshgrid(rg_sps, ld_sps);
[rggr_d2, yvgr_d2] = meshgrid(rg_sps, yv_sps);
[yvgr_d3, ldgr_d3] = meshgrid(yv_sps, ld_sps);
ld_interp_d1 = zeros(steps, steps);
ld_interp_d3 = zeros(steps, steps);
rg_interp_d1 = zeros(steps, steps);
rg_interp_d2 = zeros(steps, steps);
yv_interp_d2 = zeros(steps, steps);
yv_interp_d3 = zeros(steps, steps);
lds = squeeze(orig_data(:, 1));
rgs = squeeze(orig_data(:, 2));
yvs = squeeze(orig_data(:, 3));
vec_lds = squeeze(vecs(:, 1));
vec_rgs = squeeze(vecs(:, 2));
vec_yvs = squeeze(vecs(:, 3));
step_size_ld = abs(ld_sps(2) - ld_sps(1));
step_size_rg = abs(rg_sps(2) - rg_sps(1));
step_size_yv = abs(yv_sps(2) - yv_sps(1));

for ldc = 1:length(ld_sps)
    ld = ld_sps(ldc) + step_size_ld/2;
    for rgc = 1:length(rg_sps)
        rg = rg_sps(rgc) + step_size_rg/2;
        d1_idxs = ((lds > ld - step_size_ld/2) & (lds < ld + step_size_ld/2)) & ((rgs > rg - step_size_rg/2) & (rgs < rg + step_size_rg/2));
        
        if ~isempty(find(d1_idxs))
            ld_interp_d1(ldc, rgc) = mean(vec_lds(d1_idxs));
            rg_interp_d1(ldc, rgc) = mean(vec_rgs(d1_idxs));
        end
    end
end

for yvc = 1:length(yv_sps)
    yv = yv_sps(yvc) + step_size_yv/2;
    for rgc = 1:length(rg_sps)
        rg = rg_sps(rgc) + step_size_rg/2;
        d2_idxs = ((rgs > rg - step_size_rg/2) & (rgs < rg + step_size_rg/2)) & ((yvs > yv - step_size_yv/2) & (yvs < yv + step_size_yv/2));
        
        if ~isempty(find(d2_idxs))
            yv_interp_d2(yvc, rgc) = mean(vec_yvs(d2_idxs));
            rg_interp_d2(yvc, rgc) = mean(vec_rgs(d2_idxs));
        end
    end
end

for ldc = 1:length(ld_sps)
    ld = ld_sps(ldc) + step_size_ld/2;
    for yvc = 1:length(yv_sps)
        yv = yv_sps(yvc) + step_size_yv/2;
        d3_idxs = ((lds > ld - step_size_ld/2) & (lds < ld + step_size_ld/2)) & ((yvs > yv - step_size_yv/2) & (yvs < yv + step_size_yv/2));
        
        if ~isempty(find(d3_idxs))
            ld_interp_d3(ldc, yvc) = mean(vec_lds(d3_idxs));
            yv_interp_d3(ldc, yvc) = mean(vec_yvs(d3_idxs));
        end
    end
end

return;
end
