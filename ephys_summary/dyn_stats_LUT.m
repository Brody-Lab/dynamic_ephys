function [stats_strs, stats_aligns, stats_conds] = dyn_stats_LUT
% function [stats_strs, stats_alignments, stats_conditions, stats_includes] = stats_LUT
%
% This function builds the LUT (lookup table) for statistical comparisons used by cell_packager
%
% INPUT:
%       NONE
%
%
% OUTPUT:
%
%   stats_strs:     A list of strings that can be used to specify a
%                       stat comparison
%   stats_aligns:   integer specifying which alignment to choose from 
%                       dyn_align_LUT
%   stats_conds:    string specifying what data to use to assign conditions
%                       for each trial


align_strs = dyn_align_LUT;

%% comparisons
stats_strs{1} = 'cpokeend-choice-all';
stats_aligns(1) = strmatch('cpokeend', align_strs, 'exact');
stats_conds{1} = 'vec_data.pokedR';

stats_strs{2} = 'cpokeout-choice-all';
stats_aligns(2) = strmatch('cpokeout', align_strs, 'exact');
stats_conds{2} = 'vec_data.pokedR';

stats_strs{3} = 'cpokeout-outcome-all';
stats_aligns(3) = strmatch('cpokeout', align_strs, 'exact');
stats_conds{3} = 'vec_data.pokedR == (vec_data.bup_diff+(-.1+rand(size(vec_data.bup_diff))*.2)>0)';


