function sp = get_switch_params(which_switch, switch_params)

 sp = struct('which_switch',which_switch,...
            'clear_bad_strengths', 1, 'bad_strength', 0, 'fit_line', 1,...
            't_buffers',[0 0], 'min_pre_dur', 0, 'min_post_dur', 0,...
            'min_switch_t',0,'max_switch_t',Inf,...
            'exclude_final',0,'final_only',0,'model_smooth_wdw',100);
        
    if ~isempty(switch_params)
        %sp = p.switch_params;
        if isfield(switch_params,'which_switch')
            assert(strcmp(switch_params.which_switch, which_switch));
        end
        spfields = fieldnames(switch_params);
        for pp = 1:length(spfields)
            spfield = spfields{pp};
            sp.(spfield) = switch_params.(spfield);
        end
    end
