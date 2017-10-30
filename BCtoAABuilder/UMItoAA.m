%George Emanuel
%emanuega0@gmail.com
%#Copyright Presidents and Fellows of Harvard College, 2017.


experiment = %Sequencing experiment name
library = %Sequencing experiment library

%%

javaclasspath(%path folder containing to AdjacencyMatrix.class');

dataRoot = %root directory for data
saveRoot = %directory for output
savePath = [saveRoot library '\'];
mkdir(saveRoot);
mkdir(savePath);

%%Load fastq sequencing reads

read1File = [dataRoot library '_L001_R1_001.fastq'];
read2File = [dataRoot library '_L001_R2_001.fastq'];

read1 = fastqread(read1File);
read2 =  fastqread(read2File);

%% Extract UMIs

umi5prime = 'GCTAAGGATACGAGGCCCT';
reads = [];
goodIndices = ones(length(read1), 1);

for i=1:length(read1)
    currentSequence = read1(i).Sequence;
    currentQuality = cellfun(@(x) double(x)-64, {read1(i).Quality}, 'UniformOutput', false);
    currentQuality = currentQuality{1};
    currentSequence2 = read2(i).Sequence;
    currentQuality2 = cellfun(@(x) double(x)-64, {read2(i).Quality}, 'UniformOutput', false); 
    currentQuality2 = currentQuality2{1};
    
    umiStart = length(umi5prime) + strfind(currentSequence, umi5prime);
    reads(i).sequenceIndex = i;
    if (length(umiStart) > 0)
        reads(i).umi = currentSequence(umiStart:umiStart+19);
        reads(i).umi3prime = currentSequence(umiStart+20:umiStart+39);
        reads(i).umiStart = umiStart;
        
        reads(i).unag1 = currentSequence(umiStart + 35:end);    
        reads(i).unag1Q = currentQuality(umiStart + 35:end);
        reads(i).unag2 = currentSequence2(1:end);
        reads(i).unag2Q = currentQuality2(1:end);
    else 
        reads(i).umiStart = -1;
        goodIndices(i) = 0;
    end
end

goodReads = reads(logical(goodIndices));

goodUMIs = zeros(length(goodReads), 20);
for i=1:length(goodUMIs)
    goodUMIs(i,:) = goodReads(i).umi;
end

%% Create UMI Adjacency Matrix

%only look at unique UMIs
[uniqueUMI, ~, umiUI] = unique(goodUMIs, 'rows');


edges = AdjacencyMatrix.calculateEdges(uniqueUMI, 1) + 1;
edges = [edges; repmat(1:length(uniqueUMI), [2 1])'];


%% Find connected components of the UMI graph to find consensus UMIs, along with their consensus UnaG sequences

consensusPath = [savePath 'ConsensusSequences.mat'];

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
        currentUnaG = cell(length(currentCC), 1);
        currentUnaG2 = cell(length(currentCC), 1);
        currentQ2 = cell(length(currentCC), 1);
        
        minLength1 = 10000;
        minLength2 = 10000;
        for j=1:length(currentCC)
            currentRead = goodReads(currentCC(j));

            currentUMIs{j} = currentRead.umi;
            currentUnaG{j} = currentRead.unag1;
            currentUnaG2{j} = currentRead.unag2;
            minLength1 = min(minLength1, length(currentRead.unag1));
            minLength2 = min(minLength2, length(currentRead.unag2));
        end
        
        for j=1:length(currentCC)
            currentUnaG{j} = currentUnaG{j}(1:minLength1);
            currentUnaG2{j} = currentUnaG2{j}(1:minLength2);
        end

        %unagAlign = multialign(currentUnaG(1:min(end, 40)), 'GapOpen', 100);
        %unag2Align = multialign(currentUnaG2(1:min(end, 40)), 'GapOpen', 100);

        consensusUMIs(i).sequence = seqconsensus(currentUMIs);
        consensusUMIs(i).umiProf = seqprofile(currentUMIs, 'Alphabet', 'NT');
        consensusUMIs(i).count = length(currentCC);
        consensusUMIs(i).unag1 = currentUnaG;
        consensusUMIs(i).unag2 = currentUnaG2;
        consensusUMIs(i).unag1Prof = seqprofile(currentUnaG, 'Alphabet', 'NT');
        consensusUMIs(i).unag2Prof = seqprofile(currentUnaG2, 'Alphabet', 'NT');
    end
    
    save(consensusPath, 'consensusUMIs');
else
    load(consensusPath);
end

%% Create UMI abundance plot

abundances = [];
for i=1:length(consensusUMIs)
    abundances(i) = length(consensusUMIs(i).unag1);
end


figHandle = figure('Name', ['UMI Abundances'], 'Position', [100, 100, 250, 250]);
semilogy(abundances)
ylabel('Read count')
xlabel('UMI index')
title('Reads per UMI')
SaveFigure(figHandle, 'overwrite', true, 'formats', {'png'});


%% Extract sequencing alignment. Modify for each genetic variant studied
yfastSeq1 = 'GCAATTCAGTTGGATGGCGATGGTAACATCTTGCAGTACAATGCCGCCGAGGGTGATATTACAGGACGTGATCCCAAACAAGTGATTGGAAAAAATTTTTTCAAAGATGTAGCGCCTGGCACTGACTCACCCGAGTTTTACGGTAAGTTCAAAGAAGGCGTGGCTTCCGGTAATCTTAATACGATGTTTGAGTGGATGATCCCCACTAGCCGTGGA';
yfast3prime = 'CTTTCACCTTGGTGGGTCC';
yfast5prime = 'GGATGGTCTTGCCTTTGGAGCA';
wtAA = nt2aa(yfastSeq1, 'ALTERNATIVESTARTCODONS', 0);

for i=1:length(consensusUMIs)
    
    currentNTSequences1 = consensusUMIs(i).unag1;
    currentNTSequences2 = consensusUMIs(i).unag2;
    currentAASequences = [];
    
    for j=1:length(currentNTSequences1)
        currentAlignment1 = localalign(currentNTSequences1{j}, yfast3prime);
        currentAlignment2 = localalign(currentNTSequences2{j}, yfast5prime);

        sequence1 = '';
        if (currentAlignment1.Stop(1) - currentAlignment1.Start(1) == 18 && currentAlignment1.Start(1) > 50 && currentAlignment1.Start(1) < 70)
            sequence1 = currentNTSequences1{j}(currentAlignment1.Stop(1)-2:end);
            if length(sequence1) >= 60
                sequence1 = seqrcomplement(sequence1(1:60));
            else
                sequence1 = '';
            end
        end
        
        sequence2 = '';
        if (currentAlignment2.Stop(1) - currentAlignment2.Start(1) == 21 && currentAlignment2.Start(1) < 15)
            sequence2 = currentNTSequences2{j}(currentAlignment2.Stop(1)-2:end);
            if length(sequence2) >= 156
                sequence2 = sequence2(1:156);
            else
                sequence2 = '';
            end
        end
        
        if length(sequence1) > 0 && length(sequence2) > 0
            fullSequence = strcat(sequence2, sequence1);
            currentAASequences(end+1,1:length(wtAA)) = nt2aa(fullSequence, 'ALTERNATIVESTARTCODONS', 0, 'ACGTOnly', false);
        end
        
    end
    
    if (length(currentAASequences) > 0)

        [uA, ~, uIdx] = unique(currentAASequences, 'rows');
        [idxHist, idxI] = sort(hist(uIdx, 1:max(uIdx)));
        modeIdx = mode(uIdx);
        modeAA = uA(modeIdx, :);


        [s, a] = nwalign(wtAA, char(modeAA));

        consensusUMIs(i).aaSequence = modeAA;
        consensusUMIs(i).aaAlignment = a;
        if (length(idxHist) > 1)
            consensusUMIs(i).confidence = idxHist(end)/idxHist(end-1);
        end
    end
    i
end



%% Save UMI - aaChange table

umiAAOutput = {};
for i=1:length(consensusUMIs)
    umiAAOutput(i).umi = consensusUMIs(i).sequence;
    umiAAOutput(i).aaAlignment = consensusUMIs(i).aaAlignment;
    umiAAOutput(i).aaSequence = consensusUMIs(i).aaSequence;
    umiAAOutput(i).confidence = consensusUMIs(i).confidence;
end

umiAAPath = [savePath 'UMItoAA.mat'];
save(umiAAPath, 'umiAAOutput');

%% Confidence distribution
figHandle = figure('Name', ['confidence distribution'], 'Position', [100, 100, 250, 250]);

hist(log([umiAAOutput(:).confidence]),20)
ylabel('Counts')
xlabel('log(confidence ratio)')
title('Confidence distribution')

SaveFigure(figHandle, 'overwrite', true, 'formats', {'png'});
