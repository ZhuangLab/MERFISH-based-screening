%George Emanuel
%#Copyright Presidents and Fellows of Harvard College, 2017.

function [ bits, finalScore ] = sequenceToBitsAlign( read1, read2, bitSequences, bitPairs, spacerLength )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 5
        spacerLength = 0
    end

    

    scores1 = -ones(length(bitPairs), 1);
    scores1(:) = NaN;
    scores2 = -ones(length(bitPairs), 1);
    scores2(:) = NaN;

    for j=1:length(bitPairs)
        bitStart5 = (j-1)*(20+spacerLength) + 1;
        bitEnd5 = bitStart5 + 19;
        bcLength = (20+spacerLength)*length(bitPairs)-spacerLength;

        if (bitStart5 < length(read1))

            bitOnScore = swalign(read1, bitSequences(bitPairs(j,1)).Sequence);
            bitOffScore = swalign(read1, bitSequences(bitPairs(j,2)).Sequence); 

            scores1(j) = bitOnScore - bitOffScore;
        end

        bitStart3 = bcLength - bitEnd5 + 1;
        bitEnd3 = bitStart3 + 19;
        

        if (bitStart3 < length(read2))
            
            bitOnScore = swalign(read2, seqrcomplement(bitSequences(bitPairs(j,1)).Sequence));
            bitOffScore = swalign(read2, seqrcomplement(bitSequences(bitPairs(j,2)).Sequence)); 

            if isnan(scores1(j))
            
                scores1(j) = bitOnScore - bitOffScore;
            else
                scores1(j) = (scores1(j) + (bitOnScore - bitOffScore))*0.5;
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

