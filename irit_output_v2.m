function pc=irit_output_v2(G,DD,lab,sheet)

L=size(DD,1);
p=nan(L,1);
z1=zeros(L,6);
for k=1:L
    t=0;
    st=0;
    [p(k),t,st] = anova1(DD(k,:),lab,'off');
    c = multcompare(st,'display','off');
    [m(k,:) ss(k,:)]=grpstats(DD(k,:),lab',{'mean','sem'});    
    
    for l=1:6
        if   c(l,3)*c(l,5)>0
            z1(k,l)=sign( m(k,round(c(l,2)) ) - m(k,round(c(l,1)) ));
        end
    end
end

pc=mafdr(p,'BHFDR',1);


header={'Gene' 'anova_FDR'	'av_M+IgG' 'av_L+IgG'	'av_L+a-T' 'av_L+ a-10'...
         'std_M+IgG' 'std_L+IgG' 'std_L+a-T' 'std_L+ a-10' strcat(num2str(round(c(1,2))),'_vs_',num2str(round(c(1,1))))...
         strcat(num2str(round(c(2,2))),'_vs_',num2str(round(c(2,1)))) strcat(num2str(round(c(3,2))),'_vs_',num2str(round(c(3,1))))...
         strcat(num2str(round(c(4,2))),'_vs_',num2str(round(c(4,1)))) strcat(num2str(round(c(5,2))),'_vs_',num2str(round(c(5,1))))...
         strcat(num2str(round(c(6,2))),'_vs_',num2str(round(c(6,1))))};
                  
size(header)         
results=[G(pc<0.05) num2cell([pc(pc<0.05) m(pc<0.05,:) ss(pc<0.05,:) z1(pc<0.05,:)])];
xlswrite('test_v2',[header;results],sheet);
