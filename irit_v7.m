function irit_v7

%read and allocate data
[data,names]=xlsread('Loops Affy Data.xlsx');
D1=log2(data(:,2:65));
D_lin=data(:,2:65);
N=names(2,2:65);
G1=names(3:end,67);
[g1 g2]=grp2idx(N);
%%%%%%%%%%%%%%%%%%%%%%%


%compute stats for logged and linear data
m1=0;
for k=1:size(D1,1)
    [m1(k,1:11),s1(k,1:11)]=grpstats(D1(k,:),g1,{'mean','sem'});    
    [m1_lin(k,1:11),s1_lin(k,1:11)]=grpstats(D_lin(k,:),g1,{'mean','sem'});
end
%%%%%%%%%%%%%%%%%%%%%%%


%filtering probes with low variance
mask = genevarfilter(D1,'PERCENTILE',40);
D=D1(mask,:);
G=G1(mask);
m=m1(mask,:);
s=s1(mask,:);
m_lin=m1_lin(mask,:);
s_lin=s1_lin(mask,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%analysis at 4 hours
D_4=D(:,[find(g1==2)' find(g1==3)' find(g1==4)' find(g1==5)']);
lab_4=[g2(2*ones(1,sum((g1==2))))'  g2(3*ones(1,sum((g1==3))))' g2(4*ones(1,sum((g1==4))))' g2(5*ones(1,sum((g1==5))))'];
pc_4=irit_output_v2(G,D_4,lab_4,'4h');

%8analysis at 8 hours
D_8=D(:,[find(g1==7)' find(g1==8)' find(g1==9)' find(g1==10)']);
lab_8=[g2(7*ones(1,sum((g1==7))))'  g2(8*ones(1,sum((g1==8))))' g2(9*ones(1,sum((g1==9))))' g2(10*ones(1,sum((g1==10))))'];
pc_8=irit_output_v2(G,D_8,lab_8,'8h');

%specificity analysis at 4h
G_4=G(pc_4<0.05);
D_4l=D(pc_4<0.05,[find(g1==2)' find(g1==3)' find(g1==4)' find(g1==5)']);
lab_4l=[g2(2*ones(1,sum((g1==2))))' g2(3*ones(1,sum((g1==3))))' g2(4*ones(1,sum((g1==4))))' g2(5*ones(1,sum((g1==5))))'];
crit=grp2idx(lab_4l);

%-TNF_specific
[IDX_tnf_4 Z_tnf_4]= rankfeatures(2.^D_4l,crit==3);
G_4(Z_tnf_4>4)
%-IL10_specific
[IDX_il10_4 Z_il10_4]= rankfeatures(2.^D_4l,crit==4);
G_4(Z_il10_4>4)
 


%specificity analysis at 8h
G_8=G(pc_8<0.05);
D_8l=D(pc_8<0.05,[find(g1==7)' find(g1==8)' find(g1==9)' find(g1==10)']);
lab_8l=[g2(7*ones(1,sum((g1==7))))' g2(8*ones(1,sum((g1==8))))' g2(9*ones(1,sum((g1==9))))' g2(10*ones(1,sum((g1==10))))'];
crit=grp2idx(lab_8l);

%LPS_specific
[IDX_lps_8 Z_lps_8]= rankfeatures(D_8l,crit==2);
LPS=G_8(Z_lps_8>4);

%-TNF_specific
[IDX_tnf_8 Z_tnf_8]= rankfeatures(D_8l,crit==3);
TNF=G_8(Z_tnf_8>4);

%-IL10_specific
[IDX_il10_8 Z_il10_8]= rankfeatures(D_8l,crit==4);
IL10=G_8(Z_il10_8>4);



%%%separability plots
subplot(2,2,1);
[counts1, values1] = hist(Z_med_8, 50);
bar(values1(values1<4), counts1(values1<4),'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
hold on;
bar(values1(values1>=4), counts1(values1>=4),'facecolor','k','edgecolor','k');
xlim([-0.1 11]);

subplot(2,2,2);
[counts1, values1] = hist(Z_lps_8, 50);
bar(values1(values1<4), counts1(values1<4),'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
hold on;
bar(values1(values1>=4), counts1(values1>=4),'facecolor','k','edgecolor','k');
xlim([-0.1 11]);

subplot(2,2,3);
[counts1, values1] = hist(Z_tnf_8, 50);
bar(values1(values1<4), counts1(values1<4),'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
hold on;
bar(values1(values1>=4), counts1(values1>=4),'facecolor','k','edgecolor','k');
xlim([-0.1 11]);

subplot(2,2,4);
[counts1, values1] = hist(Z_il10_8, 50);
bar(values1(values1<4), counts1(values1<4),'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
hold on;
bar(values1(values1>=4), counts1(values1>=4),'facecolor','k','edgecolor','k');
xlim([-0.1 11]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%overlaps between signatures and factors in L-R database
[A,B]=xlsread('Lig-Rec DB 5.6.2014.xlsx');
SOLUBLE=unique(B(2:end,1:2));
SOLUBLE(1)=[];
RECEPTOR=unique(B(2:end,3:5));
RECEPTOR(1)=[];

[data_8,names_8]=xlsread('test_v2.xls','8h');
LPS_up=names_8(1+find(Z_lps_8>=4 & data_8(:,3)>max(data_8(:,[2 4 5])')'));
LPS_down=names_8(1+find(Z_lps_8>=4 & data_8(:,3)<min(data_8(:,[2 4 5])')'));

TNF_up=names_8(1+find(Z_tnf_8>=4 & data_8(:,4)>max(data_8(:,[2 3 5])')'));
TNF_down=names_8(1+find(Z_tnf_8>=4 & data_8(:,4)<min(data_8(:,[2 3 5])')'));

IL10_up=names_8(1+find(Z_il10_8>=4 & data_8(:,5)>max(data_8(:,2:4)')'));
IL10_down=names_8(1+find(Z_il10_8>=4 & data_8(:,5)<min(data_8(:,2:4)')'));

%for table 1:
Q1=intersect(SOLUBLE,LPS_down)
Q2=intersect(SOLUBLE,LPS_up)
Q3=intersect(SOLUBLE,TNF_down)
Q4=intersect(SOLUBLE,TNF_up)
Q5=intersect(SOLUBLE,IL10_down)
Q6=intersect(SOLUBLE,IL10_up)

QQ1=intersect(RECEPTOR,LPS_down)
QQ2=intersect(RECEPTOR,LPS_up)
QQ3=intersect(RECEPTOR,TNF_down)
QQ4=intersect(RECEPTOR,TNF_up)
QQ5=intersect(RECEPTOR,IL10_down)
QQ6=intersect(RECEPTOR,IL10_up)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PCA fig 3 for DC.

[data_DC,names_DC]=xlsread('Table T-mod DB .xlsx','DC');
%DC=data_DC./repmat(max(data_DC),size(data_DC,1),1)
% med=mean(data_DC(1:10,:));
% lps=mean(data_DC(11:20,:));
% lpst=mean(data_DC(21:30,:));
% lps10=mean(data_DC(31:40,:));
%[COEFF, SCORE, LATENT] = princomp(log2([med;lps;lpst;lps10]));
[COEFF, SCORE, LATENT] = princomp(log10(1+data_DC));
med=mean([SCORE(1:10,1) SCORE(1:10,2)]);
lps=mean([SCORE(11:20,1) SCORE(11:20,2)]);
lpst=mean([SCORE(21:30,1) SCORE(21:30,2)]);
lps10=mean([SCORE(31:40,1) SCORE(31:40,2)]);
DAT=[med;lps;lpst;lps10];
gscatter(DAT(:,1),DAT(:,2),[1 2 3 4],[],[],25);
cumsum(LATENT)./sum(LATENT)




%could be moved to another script...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Naive T cells
[data,names]=xlsread('PCA me please !.xlsx','data_v1');
[g1T g2T]=grp2idx(names(2:end,1));

for k=1:size(data,2)
    M(k,:) = grpstats(data(:,k),g1T)
end

[coefs, score,latent] = princomp(zscore(M'));
%scatter(score(:,1),score(:,2),100,'FaceColor','k');
gscatter(score(:,1),score(:,2),[1 2 3 4],[],[],30);
hold on;
xlim([-3.2 3.2]);
ylim([-2 2]);
cumsum(latent)./sum(latent)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Memory T cells
[data,names]=xlsread('Memory T cells for PCA.xlsx','PCA');
[g1T g2T]=grp2idx(names(2:end,1));

for k=1:size(data,2)
    M(k,:) = grpstats(data(:,k),g1T)
end

[coefs, score,latent] = princomp(zscore(M'));
%scatter(score(:,1),-score(:,2),100,'FaceColor','k');
gscatter(score(:,1),-score(:,2),[1 2 3 4],[],[],30);
hold on;
xlim([-5 5]);
ylim([-2 2]);
line([0 0],[-2 2],'LineStyle','--','Color','k')
line([-5 5],[0 0],'LineStyle','--','Color','k')
cumsum(latent)./sum(latent)




%enrichment analysis

[A,B]=xlsread('pos_reg_pro_secretion.xlsx','pos_reg_cell_comm GO 0010647');
B=unique(B(:,3));

enr_TNF=1-sum(hygecdf(length(intersect(B,TNF))-1,20000, length(B),sum(Z_tnf_8>=4)))
enr_IL10=1-sum(hygecdf(length(intersect(B,IL10))-1,20000, length(B),sum(Z_il10_8>=4))) 


%panel 4D
















TARGETS={
    'PDGFA'
    'PVR'
    'EFNB1'
    'IL20RB'
    'ITGAV'
    'PDGFRA'
    'TNFRSF1B'
    'IL21R'
    'PTPN2'};

labels={'Med','LPS','LPS+\alphaTNFR','LPS+\alphaIL10R'};
for s=1:length(TARGETS)
    k=find(strcmp(G,TARGETS(s)));
    subplot(3,3,s);
    %if length(k)==1
    [m_lin(k,7:10) s_lin(k,7:10)]
    bar(m_lin(k(1),7:10),'FaceColor','k');hold on;errorb(m_lin(k(1),7:10),s_lin(k(1),7:10));
    set(gca,'XTickLabel',{ });
    %     else
    %         bar(m_lin(k(2),7:10),'FaceColor','k');hold on;errorb(m_lin(k(2),7:10),s_lin(k(2),7:10));
    %     end
    %     if s>2
    %         set(gca,'XTickLabel',labels)
    %         rotateticklabel(gca,45);
    %     end
    %     title(TARGETS(s))
end










