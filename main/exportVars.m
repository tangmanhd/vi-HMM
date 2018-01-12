function exportVars(varsamM2, varsamsupM2, reffile, start_n, chrom) 
    if size(varsamM2,1)==0
        estvar_M2=dataset([],'VarNames','CHROM');
        fullfile(file.mainDir, file.chromSizesDir,[refGenome '.txt']);
        export(estvar_M2,'file',fullfile(pwd,'variant.txt'),'Delimiter','\t','WriteVarNames',false)
    else
        outM2=updateVars(varsamM2,varsamsupM2,reffile,start_n);
        chrom_M2=dataset([],'VarNames','CHROM');
        chrom_M2.CHROM=repmat(chrom,size(outM2,1),1);
        estvar_M2=horzcat(chrom_M2,outM2);
        export(estvar_M2,'file',fullfile(pwd,'variant.txt'),'Delimiter','\t','WriteVarNames',false)
    end
end 
