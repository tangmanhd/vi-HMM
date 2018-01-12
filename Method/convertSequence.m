function sseq = convertSequence(nseq, alphabet)
% convert numeric sequence to character sequence, written by X Wu, 2015
% input:
%         nseq: vector of any length, with possible values 1/2/3/4/5
%         alphabet: vector of length 5 or 4, nucleotide base alphabet, usually 'ACGT-' or 'ACGT'
% output:
%         sseq: vector of characters, same length as nseq
% example:
%         nseq = [2, 1, 1, 3, 4, 4]; alphabet = 'ACGT';
%         convert_sequence(nseq, alphabet)

sseq = alphabet(nseq);
