function build_dataset(rats)

if nargin < 1
    p = set_dyn_path;
    % Set which rats to analyze
    %rats  = {'H037','H066', 'H084','H129','H140'};
    %rats  = {'H130', 'H140','H141', 'H176','H190','H191'};
    rats  = p.ratlist;
end
    
p       = set_dyn_path;

if isstr(rats)
    rats = {rats};
end

p.rats  = rats;

% Do rat level analysis
for ii = 1 : length(p.rats)
    [S] = load_data(p.rats{ii},p);
    save_good_data(p.rats{ii},S,p);
    disp(' ')
end


