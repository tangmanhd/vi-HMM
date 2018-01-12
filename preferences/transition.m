function  tprob = transition(dataTypes)

%   transition build transition matrix.
%   by default dataTypes is example

%   There are four states in transition matrix: Match, SNP, Deletion,
%   Insertion
%   
    if  strcmp(dataTypes, 'example')
        tprob = [0.988, 0.008,0.002, 0.002;
             0.53 , 0.45, 0.01, 0.01;
             0.70, 0.15, 0.15 , 0.0 ;
             0.70, 0.15, 0.0 , 0.15];
    elseif strcmp(dataTypes, 'real')
        tprob = [0.9838043, 0.01474720,0.0006085089,0.0008400445;
             0.9499207, 0.04640025,0.0014855172,0.0021934910;
             0.2879631, 0.01089283,0.6994015911,0.0017424552;
             0.4163771, 0.01984721,0.0040161923,0.5597594535];
    else 
        error('Error: Data type is not valid.')
    end

end
