function [align_strs, align_args] = dyn_align_LUT(max_dur)
% function [align_strs, align_args] = align_LUT
%
% This function builds the LUT (lookup table) for PSTH "alignments" (which
% includes masking)
%
% INPUT:
%       NONE
%
% OUTPUT:
%
%   align_strs:     A list of strings that can be used to specify a
%                       particular alignment
%   align_args:     strings specifying arguments for make_rate_functions
%
% ALIGNMENT STRINGS: 
%   stimstart, stimend, premove, postmove, cpokeend, cpokeout,
%   stimend-nomask, cpokeout-nomask, stimstart-cout-mask
% 

%% alignments 
if nargin < 1
    max_dur = 2;
    fprintf('using max duration 2s\n')
end
align_strs{1} = 'stimstart';
align_args{1} = {'ref_event', 'stim_start', 'pre', 0, 'post', max_dur, 'post_mask_event', 'cpoke_end'};

align_strs{2} =  'stimend';
align_args{2} = {'ref_event', 'cpoke_end', 'pre_mask_event', 'stim_start', 'pre', max_dur, 'post_mask_event', 'cpoke_end', 'post', 0};

align_strs{3} =  'premove';
align_args{3} = {'ref_event', 'cpoke_out', 'pre_mask_event', 'stim_start', 'pre', .5+max_dur, 'post_mask_event', 'cpoke_out', 'post', 0};

align_strs{4} =  'postmove';
align_args{4} = {'ref_event', 'cpoke_out', 'pre', 0, 'post_mask_event', 'spoke_in', 'post', .5+max_dur};

align_strs{5} =  'cpokeend';
align_args{5} = {'ref_event', 'cpoke_end', 'pre', max_dur, 'post', 1, 'pre_mask_event', 'stim_start', 'post_mask_event', 'spoke_in'};

align_strs{6} =  'cpokeout';
align_args{6} = {'ref_event', 'cpoke_out', 'pre', 1+max_dur, 'post', 1, 'pre_mask_event', 'stim_start', 'post_mask_event', 'spoke_in'};

align_strs{7} =  'stimend-nomask';
align_args{7} = {'ref_event', 'cpoke_end', 'pre', max_dur-.5, 'post', 1};

align_strs{8} =  'cpokeout-nomask';
align_args{8} = {'ref_event', 'cpoke_out', 'pre', max_dur+.5, 'post', 1};

align_strs{9} =  'stimstart-cout-mask';
align_args{9} = {'ref_event', 'stim_start', 'pre', 1, 'post', max_dur, 'post_mask_event', 'cpoke_out'};

align_strs{10} =  'cout100';
align_args{10} = {'ref_event', 'cpoke_out', 'pre', 0, 'post_mask_event', 'spoke_in', 'post', 0.1};

align_strs{11} =  'stimstart-cout50';
align_args{11} = {'ref_event', 'cpoke_out', 'pre', max_dur+.5, 'pre_mask_event', 'stim_start', 'post_mask_event', 'spoke_in', 'post', 0.05};

align_strs{12} =  'spokein';
align_args{12} = {'ref_event', 'spoke_in', 'pre', 1, 'pre_mask_event', 'cpoke_out', 'post', max_dur+.5};

align_strs{13} =  'cpokein';
align_args{13} = {'ref_event', 'cpoke_start', 'pre', 2, 'post', max_dur+.5};

align_strs{14} = 'stimstart2';
align_args{14} = {'ref_event', 'stim_start', 'pre', 2, 'post', max_dur+.5};

align_strs{15} = 'lastswitch-cout-mask';
align_args{15} = {'ref_event', 'last_switch', 'pre', 1, 'post', max_dur+.5, 'post_mask_event','cpoke_out'};


