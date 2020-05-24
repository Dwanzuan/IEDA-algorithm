clear all  
clc
func_number=8;
runs=30;
SearchAgents_no=50;
max_nfes = 50000;
Data1=zeros(func_number,runs);
kj=func_number*runs;
Data2=zeros(kj,3);
Data3=zeros(func_number,runs);
t=zeros(func_number,runs);
f_mean=zeros(func_number,1);
t_mean=zeros(func_number);
E_mean=zeros(func_number);
H_mean=zeros(func_number);
cg=zeros(runs,max_nfes,func_number);
k=1;
kl=zeros(30,1);
for i=1:8%func_number
        fprintf('Problem =\t %d\n',i);
        Function_name=i;
        [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
        [fobj1]=Get_Functions_details1(Function_name);
        for j=1:runs
        fprintf('run =\t %d\n',j);
         
     [Best_score,Best_pos,cg_curve]=IEDA(SearchAgents_no,max_nfes,lb,ub,dim,fobj);
     cg(j,:,i)=cg_curve(1:50000);
              
            HP=fobj1(Best_pos);
            fprintf('\n HP = %d',HP);
            Data1(i,j)=Best_score; 
            Data2(k,:)=Best_pos()';
            Data3(i,j)=HP; 
            k=k+1;              
        end
        f_mean(i)=mean(Data1(i,:),2);% average of HPapp
        t_mean(i)=mean(t(i,:),2); %average of time
        H_mean(i)=mean(Data3(i,j),2);%average of HP
        E_mean(i)=std(Data1(i,:),0,2);% SD
      
end
