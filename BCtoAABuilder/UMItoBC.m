%George Emanuel
%emanuega0@gmail.com
%#Copyright Presidents and Fellows of Harvard College, 2017.

%% experiment specific parameters

experiment = %put root sequencing directory here
library = %put sequencing library name here

%%
javaclasspath(%path to directory containing AdjacencyMatrix.java%);

dataRoot = %Directory where all sequencing data is stored
saveRoot = %Directory where all output should be saved
savePath = [saveRoot library '\'];
mkdir(saveRoot);
mkdir(savePath);

read1File = [dataRoot library '_L001_R1_001.fastq'];
read2File = [dataRoot library '_L001_R2_001.fastq'];

read1File2 = [dataRoot library2 '_L001_R1_001.fastq'];
read2File2 = [dataRoot library2 '_L001_R2_001.fastq'];

read1 = fastqread(read1File);
read2 = fastqread(read2File);

read12 =  fastqread(read1File2);
read22 = fastqread(read2File2);

%% Extract sequence components from set 1

umi5prime = 'TTGGTGGCGGAGGTTCCTAA';
bc5prime = 'ACTGCCGCTCCTTTCCTACG';
bc3prime = 'TTGCAGTGGTAGCCGAGCAG';

reads = [];

bcLength = 320;

for i=1:length(read1)
    currentSequence = read1(i).Sequence;
    currentSequence2 = read2(i).Sequence;
    
    umiAlign = localalign(currentSequence(1:60), umi5prime);
    umiStart = umiAlign.Stop(1) + 1;
    
    bc5Align = localalign(currentSequence(1:150), bc5prime);
    bc5Start = bc5Align.Stop(1) + 1;
    
    bc3Align = localalign(currentSequence2(1:60), bc3prime);
    bc3Start = bc3Align.Stop(1) + 2;
    
    reads(i).bc5Start = bc5Start;
    reads(i).bc3Start = bc3Start;
    reads(i).umiStart = umiStart;
    
    if (length(umiStart) > 0 & bc5Start > 0 & bc3Start > 0)
        reads(i).index = i;
        reads(i).umi = currentSequence(umiStart:umiStart+19);
        reads(i).umi3prime = currentSequence(umiStart+20:umiStart+39);
        reads(i).umiStart = umiStart;
        
        reads(i).BC1 = currentSequence(bc5Start:end);   
        reads(i).BC2 = currentSequence2(bc3Start:end);
        reads(i).bc3Start = bc3Start;
        
        bc1Overlap = reads(i).BC1(bcLength - length(reads(i).BC2) + 1:end);
        bc2Overlap = seqrcomplement(reads(i).BC2(bcLength - length(reads(i).BC1)+ 1:end));
        
        if (length(bc1Overlap > 0) & length(bc2Overlap > 0))
            overlapMismatch = sum(bc1Overlap - bc2Overlap ~= 0)./length(bc1Overlap);
            reads(i).overlapMismatch = overlapMismatch;
        else
            reads(i).overlapMismatch = 1;
        end
        
        
    else 
        reads(i).umiStart = -1;
    end
end

goodReads = reads([reads.umiStart] > 0);

goodUMIs = zeros(length(goodReads), 20);
for i=1:length(goodUMIs)
    goodUMIs(i,:) = goodReads(i).umi;
end

%% Summarize offsets

figHandle = figure('Name', 'Offset distributions', 'Position', [100, 100, 750, 250]);

subplot(1,3,1)
hist([reads.umiStart], 1:60)
title('UMI Start')

subplot(1,3,2)
hist([reads.bc5Start], 1:150)
title('BC 5'' Start')

subplot(1,3,3)
hist([reads.bc3Start], 1:60)
title('BC 3'' Start')

SaveFigure(figHandle, 'overwrite', true, 'formats', {'png'});

%% Extract sequence components from set 2

umi5prime = 'TTGGTGGCGGAGGTTCCTAA';
bc5prime = 'ACTGCCGCTCCTTTCCTACG';
bc3prime1 = 'GGAGTAGTTGGTTGTTAGGA';
bc3prime2 = 'AGTTGGGTATGGAGAAAGGT';

reads2 = [];

bcLength = 320;

for i=1:length(read12)
    currentSequence = read12(i).Sequence;
    currentSequence2 = read22(i).Sequence;
    
    umiAlign = localalign(currentSequence(1:60), umi5prime);
    umiStart = umiAlign.Stop(1) + 1;
    
    bc5Align = localalign(currentSequence(1:150), bc5prime);
    bc5Start = bc5Align.Stop(1) + 1;
    
    bc3Align1 = localalign(currentSequence2(1:60), bc3prime1);
    bc3Align2 = localalign(currentSequence2(1:60), bc3prime2);
    if length(bc3Align1.Start) > 0
        if bc3Align1.Score > bc3Align2.Score
            bc3Start = bc3Align1.Start(1);
        else
            bc3Start = bc3Align2.Start(1);
        end
    else
        bc3Start = -1;
    end
    
    reads2(i).bc5Start = bc5Start;
    reads2(i).bc3Start = bc3Start;
    reads2(i).umiStart = umiStart;
    
    if (length(umiStart) > 0 & bc5Start > 0 & bc3Start > 0)
        reads2(i).index = i;
        reads2(i).umi = currentSequence(umiStart:umiStart+19);
        reads2(i).umi3prime = currentSequence(umiStart+20:umiStart+39);
        reads2(i).umiStart = umiStart;
        
        reads2(i).BC1 = currentSequence(bc5Start:end);   
        reads2(i).BC2 = currentSequence2(bc3Start:end);
        reads2(i).bc3Start = bc3Start;
        
        bc1Overlap = reads2(i).BC1(bcLength - length(reads2(i).BC2) + 1:end);
        bc2Overlap = seqrcomplement(reads2(i).BC2(bcLength - length(reads2(i).BC1)+ 1:end));
        
        if (length(bc1Overlap > 0) & length(bc2Overlap > 0))
            overlapMismatch = sum(bc1Overlap - bc2Overlap ~= 0)./length(bc1Overlap);
            reads2(i).overlapMismatch = overlapMismatch;
        else
            reads2(i).overlapMismatch = 1;
        end
        
        
    else 
        reads2(i).umiStart = -1;
    end
end

goodReads2 = reads2([reads2.umiStart] > 0);

goodUMIs2 = zeros(length(goodReads2), 20);
for i=1:length(goodUMIs2)
    goodUMIs2(i,:) = goodReads2(i).umi;
end

%% Summarize offsets 2

figHandle = figure('Name', 'Offset distributions BC2', 'Position', [100, 100, 750, 250]);

subplot(1,3,1)
hist([reads2.umiStart], 1:60)
title('UMI Start')

subplot(1,3,2)
hist([reads2.bc5Start], 1:150)
title('BC 5'' Start')

subplot(1,3,3)
hist([reads2.bc3Start], 1:60)
title('BC 3'' Start')

SaveFigure(figHandle, 'overwrite', true, 'formats', {'png'});

%% Create UMI Adjacency Matrix

%only look at unique UMIs
[uniqueUMI, ~, umiUI] = unique(goodUMIs, 'rows');


edges = AdjacencyMatrix.calculateEdges(uniqueUMI, 1) + 1;
edges = [edges; repmat(1:length(uniqueUMI), [2 1])'];

[uniqueUMI2, ~, umiUI2] = unique(goodUMIs2, 'rows');

edges2 = AdjacencyMatrix.calculateEdges(uniqueUMI2, 1) + 1;
edges2 = [edges2; repmat(1:length(uniqueUMI2), [2 1])'];

matched = AdjacencyMatrix.calculateNearestSequences(uniqueUMI2, uniqueUMI);

%% Create UMI adjacency matrix

consensusPath = [savePath 'ConsensusSequencesBC.mat'];

if ~exist(consensusPath)

    ccCountThreshold = 2;

    umiGraph = graph(edges(:,1), edges(:,2));
    umiCC = conncomp(umiGraph);
    
    umiFullCC = zeros(length(goodUMIs), 1);
    for i=1:length(umiFullCC)
        umiFullCC(i) = umiCC(umiUI(i));
    end
    
    ccCounts = hist(umiFullCC, 1:max(umiCC));
    [sortCounts, sortI] = sort(ccCounts, 'descend');
    goodCC = sortI(sortCounts > ccCountThreshold);

    consensusUMIs = [];
    for i=1:length(goodCC)
        i
        currentCC = find(umiFullCC == goodCC(i));
        currentUMIs = cell(length(currentCC), 1);
        currentBC1 = cell(length(currentCC), 1);
        currentBC2 = cell(length(currentCC), 1);
        currentQ2 = cell(length(currentCC), 1);
        
        minLength1 = 10000;
        minLength2 = 10000;
        for j=1:length(currentCC)
            currentRead = goodReads(currentCC(j));

            currentUMIs{j} = currentRead.umi;
            currentBC1{j} = currentRead.BC1;
            currentBC2{j} = currentRead.BC2;
            minLength1 = min(minLength1, length(currentRead.BC1));
            minLength2 = min(minLength2, length(currentRead.BC2));
        end
        
        for j=1:length(currentCC)
            currentBC1{j} = currentBC1{j}(1:minLength1);
            currentBC2{j} = currentBC2{j}(1:minLength2);
        end

        %unagAlign = multialign(currentUnaG(1:min(end, 40)), 'GapOpen', 100);
        %unag2Align = multialign(currentUnaG2(1:min(end, 40)), 'GapOpen', 100);

        consensusUMIs(i).sequence = seqconsensus(currentUMIs);
        consensusUMIs(i).umiProf = seqprofile(currentUMIs, 'Alphabet', 'NT');
        consensusUMIs(i).count = length(currentCC);
        consensusUMIs(i).BC1 = currentBC1;
        consensusUMIs(i).BC2 = currentBC2;
        consensusUMIs(i).bc1Prof = seqprofile(currentBC1, 'Alphabet', 'NT');
        consensusUMIs(i).bc2Prof = seqprofile(currentBC2, 'Alphabet', 'NT');
    end
    
    save(consensusPath, 'consensusUMIs');
else
    load(consensusPath);
end

%%  Add BC2 reads

    ccCountThreshold = 2;

    umiGraph = graph(edges2(:,1), edges2(:,2));
    umiCC = conncomp(umiGraph);
    
    umiFullCC = zeros(length(goodUMIs2), 1);
    for i=1:length(umiFullCC)
        umiFullCC(i) = umiCC(umiUI2(i));
    end
    
    ccCounts = hist(umiFullCC, 1:max(umiCC));
    [sortCounts, sortI] = sort(ccCounts, 'descend');
    goodCC = sortI(sortCounts > ccCountThreshold);

    minDistances = [];
    for i=1:length(goodCC)
        i
        currentCC = find(umiFullCC == goodCC(i));
        currentUMIs = cell(length(currentCC), 1);
        currentBC1 = cell(length(currentCC), 1);
        currentBC2 = cell(length(currentCC), 1);
        currentQ2 = cell(length(currentCC), 1);
        
        minLength1 = 10000;
        minLength2 = 10000;
        for j=1:length(currentCC)
            currentRead = goodReads2(currentCC(j));

            currentUMIs{j} = currentRead.umi;
            currentBC1{j} = currentRead.BC1;
            currentBC2{j} = currentRead.BC2;
            minLength1 = min(minLength1, length(currentRead.BC1));
            minLength2 = min(minLength2, length(currentRead.BC2));
        end
        
        for j=1:length(currentCC)
            currentBC1{j} = currentBC1{j}(1:minLength1);
            currentBC2{j} = currentBC2{j}(1:minLength2);
        end   
        
        consensusUMI = seqconsensus(currentUMIs);
        
        umiDistances = zeros(length(consensusUMIs),1);
        for j=1:length(umiDistances)
            umiDistances(j) = sum(consensusUMIs(j).sequence - consensusUMI ~= 0);
        end
        [minDistances(i), minIndex] = min(umiDistances);
        if minDistances(i) <= 2
            consensusUMIs(minIndex).BC12 = currentBC1;
            consensusUMIs(minIndex).BC22 = currentBC2;
            consensusUMIs(minIndex).count2 = length(currentCC);
        end
    end
    
%% Create UMI abundance plot

abundances = [];
abundances2 = [];
for i=1:length(consensusUMIs)
    abundances(i) = length(consensusUMIs(i).BC1);
    if length(consensusUMIs(i).count2) > 0
        abundances2(i) = consensusUMIs(i).count2;
    end
end


figHandle = figure('Name', ['UMI Abundances'], 'Position', [100, 100, 250, 250]);
semilogy(abundances2)
hold on;
semilogy(abundances,'r')
hold off;
ylabel('Read count')
xlabel('UMI index')
title('Reads per UMI')
SaveFigure(figHandle, 'overwrite', true, 'formats', {'png'});



%% Call barcodes

bitSequences = fastaread(%fasta file containing bit readout sequences%);

BCpairs = [1, 19;
         2, 17;
         3, 18;
         4, 22;
         5, 20;
         6, 21;
         7, 25;
         8, 23;
         9, 24;
         10, 27;
         11, 31;
         12, 26;
         13, 29;
         14, 32;
         15, 28;
         16, 33;
         30, 34;
         35, 36;
         37, 38;
         39, 40;
         41, 42;
         43, 44];

BCpairs2 = BCpairs(1:14,:);
     
for i=1:length(consensusUMIs)

    counts = consensusUMIs(i).count;
    counts2 = consensusUMIs(i).count2;
    
    i
    if (counts > 0 && length(counts2) > 0 && counts2 > 0)
        
        bitCalls = zeros(min(100,length(consensusUMIs(i).BC1)), length(BCpairs));
        bitScores = zeros(size(bitCalls));

        for j=1:min(100,length(consensusUMIs(i).BC1))
            [bitCalls(j,:), bitScores(j,:)] = sequenceToBitsAlign(consensusUMIs(i).BC1{j}, consensusUMIs(i).BC2{j}(2:end), bitSequences, BCpairs, 1);
        end            
        
        bitCalls2 = zeros(min(100,length(consensusUMIs(i).BC12)), length(BCpairs2));
        bitScores2 = zeros(size(bitCalls2));
        
        for j=1:min(100,length(consensusUMIs(i).BC12))
            [bitCalls2(j,:), bitScores2(j,:)] = sequenceToBitsAlign(consensusUMIs(i).BC12{j}, consensusUMIs(i).BC22{j}, bitSequences, BCpairs2, 1);
        end
        
        consensusUMIs(i).bitCalls = bitCalls;
        consensusUMIs(i).bitScores  = bitScores;
        
        consensusUMIs(i).bitCalls2 = bitCalls2;
        consensusUMIs(i).bitScores2  = bitScores2;
    end
    
end

%% Extract good barcodes

goodBarcodes = [];
barcodeIndex = 1;

for i=1:length(consensusUMIs)
    if (length(consensusUMIs(i).bitCalls) > 1)
        
        [uA, ~, uIdx] = unique(consensusUMIs(i).bitCalls, 'rows');
        [idxHist, idxI] = sort(hist(uIdx, 1:max(uIdx)));
        modeIdx = idxI(end);
        consensusMode = uA(modeIdx, :);
        consensusMode(isnan(sum(consensusUMIs(i).bitScores))) = nan;
        
        
        [uA2, ~, uIdx2] = unique(consensusUMIs(i).bitCalls2, 'rows');
        [idxHist2, idxI2] = sort(hist(uIdx2, 1:max(uIdx2)));
        modeIdx2 = idxI2(end);
        consensusMode2 = uA2(modeIdx2, :);
        consensusMode2(numel(consensusMode2)+1:numel(consensusMode)) = NaN;
        
        
        calledBarcode = zeros(length(consensusMode), 1);
        for j=1:length(calledBarcode)
            if ~isnan(consensusMode2(j))
                calledBarcode(j) = consensusMode2(j);
            elseif ~isnan(consensusMode(j))
                calledBarcode(j) = consensusMode(j);
            end
            
            if ~isnan(consensusMode(j)) && ~isnan(consensusMode2(j))
                if consensusMode(j) ~= consensusMode2(j)
                    disp(['Barcodes do not match ' num2str(j) ', ' num2str(i)])
                end
            end
        end
        
        goodBarcodes(barcodeIndex).bc = calledBarcode;
        goodBarcodes(barcodeIndex).umi = consensusUMIs(i).sequence;
        if (length(idxHist) > 1)
            goodBarcodes(barcodeIndex).confidence = idxHist(end)/idxHist(end-1);
        end
        barcodeIndex = barcodeIndex + 1;
    end
end


goodBCPath = [savePath '\UMItoBC.mat'];
save(goodBCPath, 'goodBarcodes')


%% Confidence distribution
figHandle = figure('Name', ['confidence distribution'], 'Position', [100, 100, 250, 250]);

hist(log([goodBarcodes(:).confidence]),20)
ylabel('Counts')
xlabel('log(confidence ratio)')
title('Confidence distribution')

SaveFigure(figHandle, 'overwrite', true, 'formats', {'png'});
