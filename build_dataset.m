function build_dataset(rats)

if nargin < 1
    % Set which rats to analyze
    %rats  = {'H037','H066', 'H084','H129','H140'};
    %rats  = {'H130', 'H140','H141', 'H176','H190','H191'};
    rats  = {'H129'};
end
    
p       = set_dyn_path;
p.rats  = rats;

% Do rat level analysis
for ii = 1 : length(p.rats)
    [S] = load_data(p.rats{ii},p);
    save_good_data(p.rats{ii},S,p);
    disp(' ')
end


