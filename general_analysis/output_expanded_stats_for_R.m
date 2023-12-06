clear all; clc;

%% organize the data!

load('../data/obsstats.mat', 'data_table');
writetable(data_table, '../data/data_table.csv');

obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});
obs_data = removevars(obs_data, {'obs_name'});
pop_data = grpstats(obs_data, {'exp_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});

%% output for R

writetable(obs_data, '../data/obsstats_for_R.csv');
writetable(pop_data, '../data/popstats_for_R.csv');

%% do it again fro the dark_exclusion data
load('../data/dark_exc_obsstats.mat', 'data_table');
writetable(data_table, '../data/data_table_dark_exc.csv');

obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});
obs_data = removevars(obs_data, {'obs_name'});
pop_data = grpstats(obs_data, {'exp_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});

%% output for R

writetable(obs_data, '../data/obsstats_for_R_dark_exc.csv');
writetable(pop_data, '../data/popstats_for_R_dark_exc.csv');

%% do it one more time for the white walls mask data
load('../data/white_walls_obsstats.mat', 'data_table');
writetable(data_table, '../data/data_table_white_walls.csv');

obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});
obs_data = removevars(obs_data, {'obs_name'});
pop_data = grpstats(obs_data, {'exp_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});

%% output for R

writetable(obs_data, '../data/obsstats_for_R_white_walls.csv');
writetable(pop_data, '../data/popstats_for_R_white_walls.csv');