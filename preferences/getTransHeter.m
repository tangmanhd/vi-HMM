function [tprob, hetrate] = getTransHeter(dataTypes)
%  getTransHeter function to get transiton probability matrix and
%  heterozygous rate 


    tprob = transition(dataTypes);

    if strcmp(dataTypes,'example')
        hetrate = 0.01;
    elseif strcmp(dataTypes, 'real')
        hetrate = 0.001;
    end
end

