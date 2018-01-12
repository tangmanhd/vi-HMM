function runHMMvar(dataTypes,transMatrix, heteroRate, refname, samfilename, diploid) 

    if ~exist('dataTypes','var')
        dataTypes = 'example';
    end

    if  ~exist('transMatrix','var')
        [tprob, ~] = getTransHeter (dataTypes);
    else
        tprob = transMatrix;
    end
    
    if  ~exist('heteroRate','var')
        [~, hetrate] = getTransHeter (dataTypes);
    else
        hetrate = heteroRate;
    end
    
    if ~exist('refname','var')
        refname = 'ref.fa';
    end
    
    if ~exist('samfilename','var')
        samfilename = 'example.sam';
    end
    
    if ~exist('diploid','var')
        diploid = 2;
    end
    

ref.filePath = getrefPaths(refname);
[~, reffile] = fastaread(ref.filePath.refFile);

[dat, refseq_tru, aln_read, aln_quality, aln_refindex] = getAlignMatrixMain(samfilename, 1, reffile);

[varsam,varsamsup] = runVarCalling(aln_read, aln_quality, dat, aln_refindex, refseq_tru, diploid, [], tprob, hetrate);

exportVars(varsam, varsamsup, reffile, 1, 'simref');
end


