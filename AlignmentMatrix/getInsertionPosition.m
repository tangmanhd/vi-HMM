function [sam_I]= findinsloc(dati)

    cig_index = regexp(dati.CigarString,'[A-Z]');
    cig_num =str2double(strsplit(regexprep(dati.CigarString,'[A-Z]',' ')));

    cig_index_I = regexp(dati.CigarString,'I');
    cig_num_I = cig_num;
    cig_index_S = regexp(dati.CigarString,'S');

    if size(cig_index_I,2)~=0
        for t=1:length(cig_index_S)
            S_loc=find(cig_index==cig_index_S(t));
            cig_num_I(S_loc) = 0; 
        end
        for t=1:length(cig_index_I)
        I_loc=find(cig_index==cig_index_I(t));
        cig_num_I(I_loc) = 0; 
        sam_I(t,1)= dati.Position+sum(cig_num_I(1:I_loc));
        sam_I(t,2)= cig_num(I_loc);
        end  
    else
        sam_I = [0,0];
    end
end

 