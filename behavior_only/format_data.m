function [data] = format_data(S)

c = 1;
for i=1:length(S.sessid)
    for j=1:S.n_done(i);

    data(c).T           = S.pd{i}.bupsdata{j}.T;
    data(c).leftbups    = S.pd{i}.bupsdata{j}.left;
    data(c).rightbups   = S.pd{i}.bupsdata{j}.right;
    data(c).Delta       =  length(data(c).rightbups) - length(data(c).leftbups);
    data(c).hit         = S.pd{i}.hits(j);
    data(c).gamma       = S.pd{i}.bupsdata{j}.gamma;
    data(c).pokedR      = (data(c).hit & S.pd{i}.bupsdata{j}.correctAnswer) | (~data(c).hit & ~S.pd{i}.bupsdata{j}.correctAnswer);
    data(c).Hazard      = S.pd{i}.bupsdata{j}.hazard;
    data(c).sessiondate = S.sessiondate{i}; 
    data(c).sessid      = S.sessid(i);
    data(c).genEndState = S.pd{i}.bupsdata{j}.genEndState;
    data(c).genSwitchTimes = S.pd{i}.bupsdata{j}.genSwitchTimes;
    data(c).correctAnswer =  S.pd{i}.bupsdata{j}.correctAnswer;
    data(c).evidenceRatio = S.pd{i}.bupsdata{j}.evidenceRatio;
    c= c + 1;

    end
end

