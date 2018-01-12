function file = getrefPaths(refname)
%   getrefPaths Build paths for reference file.
%   Input is reference genome build used. Default is ref.fa
%
%   Author: Man Tang

if nargin == 0
    refGenome = 'ref.fa';
else
    refGenome = refname;
end

file.mainDir = pwd;

file.refDir = 'data';

file.refFile = fullfile(file.mainDir, file.refDir,refGenome);

end
