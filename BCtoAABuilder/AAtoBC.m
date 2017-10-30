%George Emanuel
%emanuega0@gmail.com
%#Copyright Presidents and Fellows of Harvard College, 2017.

experiment = %sequencing experiment name
bcLibrary = %sequencing library for barcode
aaLibrary = %sequencing library for genetic variant


%%
saveRoot = %path to save analysis output
aaFile = [saveRoot aaLibrary '\UMItoAA.mat'];
bcFile = [saveRoot bcLibrary '\UMItoBC.mat'];

load(aaFile);
load(bcFile);

%%Pair barcode and genetic variant output

bcToAA = [];

for i=1:length(umiAAOutput)
    
    umiDistances = zeros(length(goodBarcodes), 1);
    searchUMI = seqrcomplement(umiAAOutput(i).umi);
    for j=1:length(umiDistances)
        umiDistances(j) = sum(goodBarcodes(j).umi - searchUMI ~= 0);
    end
    
    [minDistance, bcIndex] = min(umiDistances);
   
    if minDistance < 2
        
        bcToAA(end+1).umi = umiAAOutput(i).umi;
        bcToAA(end).aaSequence = umiAAOutput(i).aaSequence;
        bcToAA(end).aaAlignment = umiAAOutput(i).aaAlignment;
        bcToAA(end).barcode = goodBarcodes(bcIndex).bc';
        
        if length(umiAAOutput(i).confidence) > 0
            bcToAA(end).aaConfidence = umiAAOutput(i).confidence;
        else
            bcToAA(end).aaConfidence = 0;
        end
        
        
        if length(goodBarcodes(bcIndex).confidence) > 0
            bcToAA(end).bcConfidence = goodBarcodes(bcIndex).confidence;
        else
            bcToAA(end).bcConfidence = 0;
        end
    end
end

bcToAAFile = [saveRoot aaLibrary '\BCtoAA.mat'];

save(bcToAAFile, 'bcToAA');



