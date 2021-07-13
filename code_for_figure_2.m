function code_for_figure_2()

%%%%%%%%%%%
%FIGURE 2E%
%%%%%%%%%%%

load('tetrodeRecordings_OC_2s.mat')
figure
norm_type = 'probability';
lim_hist = [0:1:15];
lim_nb_od = -1;

%Number of odors eliciting ipsi response (Figure 2E)
M = [C{1};C{2};C{3}];
Mi = M(:,1:15,2);
M = Mi;
M = sum(M,2);
M = M(find(M~=lim_nb_od));
subplot(1,2,1)
H = histogram(M,lim_hist,'Normalization',norm_type,'DisplayStyle','stairs','EdgeColor','r');

M = [C{4};C{5};C{6};C{7}];
Mi = M(:,1:15,2);
M = Mi;
M = sum(M,2);
M = M(find(M~=lim_nb_od));
subplot(1,2,2)
H = histogram(M,lim_hist,'Normalization',norm_type,'DisplayStyle','stairs','EdgeColor','b');

%Number of odors eliciting contra response (Figure 2E)
M = [C{1};C{2};C{3}];
Mc = M(:,1:15,1);
M = Mc;
M = sum(M,2);
M = M(find(M~=lim_nb_od));
subplot(1,2,1)
hold on
histogram(M,lim_hist,'Normalization',norm_type,'DisplayStyle','stairs','EdgeColor','m');

M = [C{4};C{5};C{6};C{7}];
Mc = M(:,1:15,1);
M = Mc;
M = sum(M,2);
M = M(find(M~=lim_nb_od));
subplot(1,2,2)
hold on
histogram(M,lim_hist,'Normalization',norm_type,'DisplayStyle','stairs','EdgeColor','c');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 2G%
%%%%%%%%%%%

%Positive-nagative balance (Figure 2G)
load('tetrodeRecordings_OC_2s.mat')
Results = {};
for mouse = 1:7
    M = mean(A{mouse},4);
    mm = zeros(size(M,1),15);
    Msign = C{mouse};
    for neuron = 1:size(M,1)
        for od = 1:15
            mm(neuron,od,1) = (M(neuron,od,2)- M(neuron,16,2))*Msign(neuron,od,2);
            mm(neuron,od,2) = (M(neuron,od,1)- M(neuron,16,1))*Msign(neuron,od,1);
        end
    end
    Results{mouse} = mm;
end
Raoni = [];
Rapci = [];
Rppci = [];
for mouse = 1:7
    M = Results{mouse};
    Mr = M(:,:,1);
    Ml = M(:,:,2);
    for n = 1:size(M,1)
        m = [Mr(n,:)];
        pos = length(find(m~=0 & m>0));
        neg = length(find(m~=0 & m<0));
        if pos~=0 || neg~=0
            if mouse<=3
                Raoni = [Raoni,100*pos/(pos+neg)];
            else
                Rapci = [Rapci,100*pos/(pos+neg)];
            end
        end
    end
end
Raonc = [];
Rapcc = [];
for mouse = 1:7
    M = Results{mouse};
    Mr = M(:,:,1);
    Ml = M(:,:,2);
    for n = 1:size(M,1)
        m = [Ml(n,:)];
        pos = length(find(m~=0 & m>0));
        neg = length(find(m~=0 & m<0));
        if pos~=0 || neg~=0
            if mouse<=3
                Raonc = [Raonc,100*pos/(pos+neg)];
            else
                Rapcc = [Rapcc,100*pos/(pos+neg)];
            end
        end
    end
end

figure
hlim = [0:20:100];
nm = 'probability';
subplot(1,2,1)
hold on
histogram(Raoni,hlim,'Normalization',nm,'DisplayStyle','stairs','EdgeColor','r')
subplot(1,2,2)
hold on
histogram(Rapci,hlim,'Normalization',nm,'DisplayStyle','stairs','EdgeColor','b')
subplot(1,2,1)
hold on
histogram(Raonc,hlim,'Normalization',nm,'DisplayStyle','stairs','EdgeColor','m')
subplot(1,2,2)
hold on
histogram(Rapcc,hlim,'Normalization',nm,'DisplayStyle','stairs','EdgeColor','c')
ylim([0 0.6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 2F%
%%%%%%%%%%%

load('tetrodeRecordings_OC_2s.mat')
figure
hlim = 0:1:25;
nor = 'probability';

%Response magnitude to ipsi presentations (Figure 2F)
val = [];
cat = [];
for i = 1:7
    X = [];
    for odor_id = 1:15
        for side = 2
            Mouse = A{i};
            WhichSignif = C{i};
            WhichSignif = squeeze(WhichSignif(:,odor_id,side));
            WhichSignif = find(WhichSignif>=1);
            Blank = mean(squeeze(Mouse(WhichSignif,16,side,:)),2);
            AllTrials = mean(squeeze(Mouse(WhichSignif,odor_id,side,:)),3);
            for n = 1:length(WhichSignif)
                AllTrials(n,:) = AllTrials(n,:) - Blank(n);
            end
            AllTrials = mean(AllTrials,2);
            X = [X;AllTrials];
        end
    end
    val = [val;abs(X)];
    if i<=3
        cat = [cat;1*ones(length(X),1)];
    else
        cat = [cat;2*ones(length(X),1)];
    end
end
subplot(1,2,1)
hold on
histogram(val(find(cat==1)),hlim,'Normalization',nor,'DisplayStyle','stairs','EdgeColor','r')
xlim([0 10])
subplot(1,2,2)
hold on
histogram(val(find(cat==2)),hlim,'Normalization',nor,'DisplayStyle','stairs','EdgeColor','b')
xlim([0 10])

%Response magnitude to contra presentations (Figure 2F)
val2 = [];
cat2 = [];
for i = 1:7
    X = [];
    for odor_id = 1:15
        for side = 1
            Mouse = A{i};
            WhichSignif = C{i};
            WhichSignif = squeeze(WhichSignif(:,odor_id,side));
            WhichSignif = find(WhichSignif>=1);
            Blank = mean(squeeze(Mouse(WhichSignif,16,side,:)),2);
            AllTrials = mean(squeeze(Mouse(WhichSignif,odor_id,side,:)),3);
            for n = 1:length(WhichSignif)
                AllTrials(n,:) = AllTrials(n,:) - Blank(n);
            end
            AllTrials = mean(AllTrials,2);
            X = [X;AllTrials];
        end
    end
    val2 = [val2;abs(X)];
    if i<=3
        cat2 = [cat2;1*ones(length(X),1)];
    else
        cat2 = [cat2;2*ones(length(X),1)];
    end
end
subplot(1,2,1)
hold on
histogram(val2(find(cat2==1)),hlim,'Normalization',nor,'DisplayStyle','stairs','EdgeColor','m')
xlim([0 10])
subplot(1,2,2)
hold on
histogram(val2(find(cat2==2)),hlim,'Normalization',nor,'DisplayStyle','stairs','EdgeColor','c')
xlim([0 10])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 2B%
%%%%%%%%%%%

%Basal firing rate, per region (Figure 2B)
load('tetrodeRecordings_OC_2s.mat')
whichRegion = [1 1 1 2 2 2 2];
allBasalFR = [];
allCat = [];
for region = 1:2
    for mouse = find(whichRegion==region)
        allBasalFR = [allBasalFR ; BasalFR{mouse}];
        allCat = [allCat ; region*ones(length(A{mouse}),1)];
    end
end
figure
boxplot(allBasalFR,allCat,'PlotStyle','compact')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 2D%
%%%%%%%%%%%

load('tetrodeRecordings_OC_2s.mat')

%Response matrix, example from mouse AON1 (Figure 2D left)
mouse = 1;
Examples = [10 37 36 1 70 106 29 17 6 23]; %neurons picked randomly
for odor_id = 1:15
    for mouse_id = mouse
        Mouse = A{mouse_id};
        WhichSignif = C{mouse_id};
        WhichSignif = squeeze(WhichSignif(:,odor_id,:));
        neuron_id = Examples;
        WhichSignif = WhichSignif(neuron_id,:);
        Blank_R = mean(squeeze(Mouse(:,16,2,:)),2);
        Blank_L = mean(squeeze(Mouse(:,16,1,:)),2);
        Rall = squeeze(Mouse(neuron_id,odor_id,2,:));
        Lall = squeeze(Mouse(neuron_id,odor_id,1,:));
        WhichSignif = ones(size(WhichSignif));
        for n = 1:length(neuron_id)
            Rall(n,:) = (Rall(n,:) - Blank_R(neuron_id(n)))*WhichSignif(n,2);
            Lall(n,:) = (Lall(n,:) - Blank_L(neuron_id(n)))*WhichSignif(n,1);
        end
    end
end
Results = {};
for mouse_id = mouse
    RRes = [];
    LRes = [];
    BRes = [];
    RSEM = [];
    LSEM = [];
    for odor_id = 1:15
        Mouse = A{mouse_id};
        WhichSignif = C{mouse_id};
        WhichSignif = squeeze(WhichSignif(:,odor_id,:));
        WhichSignif = WhichSignif|1; %%%comment this line to display significant responses only
        Rres = zeros([length(Mouse),1]);
        Lres = zeros([length(Mouse),1]);
        Bres = zeros([length(Mouse),1]);
        Rsem = zeros([length(Mouse),1]);
        Lsem = zeros([length(Mouse),1]);
        Blank_R = squeeze(Mouse(:,16,2,:));
        Blank_L = squeeze(Mouse(:,16,1,:));
        Rall = squeeze(Mouse(:,odor_id,2,:));
        Lall = squeeze(Mouse(:,odor_id,1,:));
        for n = 1:length(Mouse)
            Rres(n) = WhichSignif(n,2)*(mean(Rall(n,:)) - mean(Blank_R(n,:)));
            Lres(n) = WhichSignif(n,1)*(mean(Lall(n,:)) - mean(Blank_L(n,:)));
            rrr = Rall(n,:) - mean(Blank_R(n,:));
            lll = Lall(n,:) - mean(Blank_L(n,:));
            Rsem(n) = std(rrr)/sqrt(length(rrr));
            Lsem(n) = std(lll)/sqrt(length(lll));
        end
        RRes = [RRes,Rres(Examples)];
        LRes = [LRes,Lres(Examples)];
        RSEM = [RSEM,Rsem(Examples)];
        LSEM = [LSEM,Rsem(Examples)];
    end
    Results{mouse_id,1} = RRes;
    Results{mouse_id,2} = LRes;
    Results{mouse_id,3} = RSEM;
    Results{mouse_id,4} = LSEM;
end

Ipsi = Results{mouse_id,1};
Contra = Results{mouse_id,2};
IpsiSEM = Results{mouse_id,3};
ContraSEM = Results{mouse_id,4};

[Ipsi,I] = sort(Ipsi,2);
for neuron = 1:10
    Contra(neuron,:) = Contra(neuron,I(neuron,:));
    IpsiSEM(neuron,:) = IpsiSEM(neuron,I(neuron,:));
    ContraSEM(neuron,:) = ContraSEM(neuron,I(neuron,:));
end
Contra = Contra(:,15:-1:1);
ContraSEM = ContraSEM(:,15:-1:1);

M = [Ipsi Contra]';
Msem = [IpsiSEM ContraSEM]';
isBilCorr = [1 2 5 6 7];
figure
for n = 1:10
    subplot(5,2,n)
    hold on
    toPlot = M(:,n);
    toPlotError = Msem(:,n);
    if isempty(find(isBilCorr==n))
        errorbar([-15:-1 1:15],toPlot,toPlotError,'k-')
        plot([-15:-1 1:15],toPlot,'ro')
    else
        errorbar([-15:-1 1:15],toPlot,toPlotError,'b-')
        plot([-15:-1 1:15],toPlot,'ro')
    end
    yLim = [ceil(min(toPlot))-1,ceil(max(toPlot))];
    plot([0 0],yLim,'--r')
    grid on
%     ylim(yLim)
    xlim([-16 16])
end

whichSignif = C{mouse};
whichSignif = whichSignif(Examples,1:15,:);
whichSignifIpsi = whichSignif(:,:,1);
whichSignifContra = whichSignif(:,:,2);
for neuron = 1:10;
    whichSignifIpsi(neuron,:) = whichSignifIpsi(neuron,I(neuron,:));
    whichSignifContra(neuron,:) = whichSignifContra(neuron,I(neuron,:));
end
whichSignifContra = whichSignifContra(:,15:-1:1);
M = [whichSignifIpsi whichSignifContra]';
figure
for n = 1:10
    subplot(5,2,n)
    hold on
    toPlot = M(:,n);
    if isempty(find(isBilCorr==n))
        plot([-15:-1 1:15],toPlot,'k-o')
    else
        plot([-15:-1 1:15],toPlot,'b-o')
    end
    yLim = [ceil(min(toPlot))-1,ceil(max(toPlot))];
    plot([0 0],yLim,'--r')
    grid on
    ylim(yLim)
    xlim([-16 16])
end

%Response matrix, example from mouse APC1 (Figure 2D right)
mouse = 4;
Examples = [222 246 35 140 170 27 75 145 264 307]; %neurons picked randomly
for odor_id = 1:15
    for mouse_id = mouse
        Mouse = A{mouse_id};
        WhichSignif = C{mouse_id};
        WhichSignif = squeeze(WhichSignif(:,odor_id,:));
        neuron_id = Examples;
        WhichSignif = WhichSignif(neuron_id,:);
        Blank_R = mean(squeeze(Mouse(:,16,2,:)),2);
        Blank_L = mean(squeeze(Mouse(:,16,1,:)),2);
        Rall = squeeze(Mouse(neuron_id,odor_id,2,:));
        Lall = squeeze(Mouse(neuron_id,odor_id,1,:));
        WhichSignif = ones(size(WhichSignif));
        for n = 1:length(neuron_id)
            Rall(n,:) = (Rall(n,:) - Blank_R(neuron_id(n)))*WhichSignif(n,2);
            Lall(n,:) = (Lall(n,:) - Blank_L(neuron_id(n)))*WhichSignif(n,1);
        end
    end
end
Results = {};
for mouse_id = mouse
    RRes = [];
    LRes = [];
    BRes = [];
    RSEM = [];
    LSEM = [];
    for odor_id = 1:15
        Mouse = A{mouse_id};
        WhichSignif = C{mouse_id};
        WhichSignif = squeeze(WhichSignif(:,odor_id,:));
        WhichSignif = WhichSignif|1; %%%comment this line to display significant responses only
        Rres = zeros([length(Mouse),1]);
        Lres = zeros([length(Mouse),1]);
        Bres = zeros([length(Mouse),1]);
        Rsem = zeros([length(Mouse),1]);
        Lsem = zeros([length(Mouse),1]);
        Blank_R = squeeze(Mouse(:,16,2,:));
        Blank_L = squeeze(Mouse(:,16,1,:));
        Rall = squeeze(Mouse(:,odor_id,2,:));
        Lall = squeeze(Mouse(:,odor_id,1,:));
        for n = 1:length(Mouse)
            Rres(n) = WhichSignif(n,2)*(mean(Rall(n,:)) - mean(Blank_R(n,:)));
            Lres(n) = WhichSignif(n,1)*(mean(Lall(n,:)) - mean(Blank_L(n,:)));
            rrr = Rall(n,:) - mean(Blank_R(n,:));
            lll = Lall(n,:) - mean(Blank_L(n,:));
            Rsem(n) = std(rrr)/sqrt(length(rrr));
            Lsem(n) = std(lll)/sqrt(length(lll));
        end
        RRes = [RRes,Rres(Examples)];
        LRes = [LRes,Lres(Examples)];
        RSEM = [RSEM,Rsem(Examples)];
        LSEM = [LSEM,Rsem(Examples)];
    end
    Results{mouse_id,1} = RRes;
    Results{mouse_id,2} = LRes;
    Results{mouse_id,3} = RSEM;
    Results{mouse_id,4} = LSEM;
end

Ipsi = Results{mouse_id,1};
Contra = Results{mouse_id,2};
IpsiSEM = Results{mouse_id,3};
ContraSEM = Results{mouse_id,4};

[Ipsi,I] = sort(Ipsi,2);
for neuron = 1:10
    Contra(neuron,:) = Contra(neuron,I(neuron,:));
    IpsiSEM(neuron,:) = IpsiSEM(neuron,I(neuron,:));
    ContraSEM(neuron,:) = ContraSEM(neuron,I(neuron,:));
end
Contra = Contra(:,15:-1:1);
ContraSEM = ContraSEM(:,15:-1:1);

M = [Ipsi Contra]';
Msem = [IpsiSEM ContraSEM]';
isBilCorr = [1 2 5 6 7];
figure
for n = 1:10
    subplot(5,2,n)
    hold on
    toPlot = M(:,n);
    toPlotError = Msem(:,n);
    if isempty(find(isBilCorr==n))
        errorbar([-15:-1 1:15],toPlot,toPlotError,'k-')
        plot([-15:-1 1:15],toPlot,'ro')
    else
        errorbar([-15:-1 1:15],toPlot,toPlotError,'b-')
        plot([-15:-1 1:15],toPlot,'ro')
    end
    yLim = [ceil(min(toPlot))-1,ceil(max(toPlot))];
    plot([0 0],yLim,'--r')
    grid on
%     ylim(yLim)
    xlim([-16 16])
end

whichSignif = C{mouse};
whichSignif = whichSignif(Examples,1:15,:);
whichSignifIpsi = whichSignif(:,:,1);
whichSignifContra = whichSignif(:,:,2);
for neuron = 1:10;
    whichSignifIpsi(neuron,:) = whichSignifIpsi(neuron,I(neuron,:));
    whichSignifContra(neuron,:) = whichSignifContra(neuron,I(neuron,:));
end
whichSignifContra = whichSignifContra(:,15:-1:1);
M = [whichSignifIpsi whichSignifContra]';
figure
for n = 1:10
    subplot(5,2,n)
    hold on
    toPlot = M(:,n);
    if isempty(find(isBilCorr==n))
        plot([-15:-1 1:15],toPlot,'k-o')
    else
        plot([-15:-1 1:15],toPlot,'b-o')
    end
    yLim = [ceil(min(toPlot))-1,ceil(max(toPlot))];
    plot([0 0],yLim,'--r')
    grid on
    ylim(yLim)
    xlim([-16 16])
end