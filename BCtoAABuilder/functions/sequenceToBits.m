%George Emanuel
%#Copyright Presidents and Fellows of Harvard College, 2017.

function [ bits, finalScore ] = sequenceToBits( read1, read2, bitSequences, bitPairs, spacerLength )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 5
        spacerLength = 0;
    end

    if (size(read1, 1) == 1)
        bc1Prof = seqprofile(read1, 'Alphabet', 'NT');
    else
        bc1Prof = read1;
    end
    
    if (size(read2, 1) == 1)
        bc2Prof = seqprofile(read2, 'Alphabet', 'NT');
    else
        bc2Prof = read2;
    end
    

    scores1 = -ones(length(bitPairs), 1);
    scores1(:) = NaN;
    scores2 = -ones(length(bitPairs), 1);
    scores2(:) = NaN;

    for j=1:length(bitPairs)
        bitStart5 = (j-1)*(20+spacerLength) + 1;
        bitEnd5 = bitStart5 + 19;
        bcLength = (20+spacerLength)*length(bitPairs)-spacerLength;

        if (bitStart5 < length(bc1Prof))
            if (bitEnd5 < length(bc1Prof))
                bit1 = bc1Prof(:,bitStart5:bitEnd5);

                bitOnProf = seqprofile(bitSequences(bitPairs(j,1)).Sequence, 'Alphabet', 'NT');
                bitOffProf = seqprofile(bitSequences(bitPairs(j,2)).Sequence, 'Alphabet', 'NT');

                scores1(j) = sum(sum(bit1.*bitOnProf)) - sum(sum(bit1.*bitOffProf));
            else
                bit1 = bc1Prof(:,bitStart5:end);
               
                bitOnProf = seqprofile(bitSequences(bitPairs(j,1)).Sequence(1:size(bit1,2)), 'Alphabet', 'NT');
                bitOffProf = seqprofile(bitSequences(bitPairs(j,2)).Sequence(1:size(bit1,2)), 'Alphabet', 'NT');
                
                scores1(j) = sum(sum(bit1.*bitOnProf)) - sum(sum(bit1.*bitOffProf));
            end
        end

        bitStart3 = bcLength - bitEnd5 + 1;
        bitEnd3 = bitStart3 + 19;
        

        if (bitStart3 < length(bc2Prof))
            if (bitEnd3 < length(bc2Prof))
                bit2 = bc2Prof(:,bitStart3:bitEnd3);

                bitOnProf = seqprofile(seqrcomplement(bitSequences(bitPairs(j,1)).Sequence), 'Alphabet', 'NT');
                bitOffProf = seqprofile(seqrcomplement(bitSequences(bitPairs(j,2)).Sequence), 'Alphabet', 'NT');


                scores2(j) = sum(sum(bit2.*bitOnProf)) - sum(sum(bit2.*bitOffProf));
            else
                bit2 = bc2Prof(:,bitStart3:end);

                bitOnProf = seqprofile(seqrcomplement(bitSequences(bitPairs(j,1)).Sequence), 'Alphabet', 'NT');
                bitOffProf = seqprofile(seqrcomplement(bitSequences(bitPairs(j,2)).Sequence), 'Alphabet', 'NT');

                bitOnProf = bitOnProf(:, 1:size(bit2,2));
                bitOffProf = bitOffProf(:, 1:size(bit2,2));
                
                
                scores2(j) = sum(sum(bit2.*bitOnProf)) - sum(sum(bit2.*bitOffProf));
            end
        end
    end

        
    finalScore = zeros(length(bitPairs), 1);
    for j=1:length(finalScore)
        sumCount = 0;
        sumValue = 0;
        if (~isnan(scores1(j)))
            sumValue = sumValue + scores1(j);
            sumCount = sumCount + 1;
        end
        if (~isnan(scores2(j)))
            sumValue = sumValue + scores2(j);
            sumCount = sumCount + 1;
        end
        finalScore(j) = sumValue/sumCount;
    end
    
    bits = finalScore > 0;
    
end

