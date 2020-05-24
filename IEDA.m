%___________________________________________________________________%
%  IEDA Algorithm  source codes demo version 1.0                    %
%                                                                   %
%  Developed in MATLAB R2017b                                       %
%                                                                   %
%  Author and programmer: Dai wanxuan et.al                         %
%                                                                   %
%         e-Mail:1187794879@qq.com                                  %
%                                                                   %
%       URL:https://github.com/Dwanzuan/IEDA-algorithm              %
%                                                                   %

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run IEDA: [Best_score,Best_pos,cg_curve]=DA(SearchAgents_no,max_nfes,lb,ub,dim,fobj)
%__________________________________________
function [Best_score,Best_pos,cg_curve]=IEDA(SearchAgents_no,max_nfes,lb,ub,dim,fobj)
display('IEDA is optimizing your problem');
rand('seed', sum(100 * clock)); %create 15 digit after decimal random number
cg_curve=zeros(1,max_nfes);
Delta_max=(ub-lb)/10;
DeltaX=zeros(dim,SearchAgents_no);
X=initialization(SearchAgents_no,dim,ub,lb);
Neighbours_DeltaX=zeros(dim,SearchAgents_no);
Neighbours_X=zeros(dim,SearchAgents_no);
Fitness=zeros(1,SearchAgents_no);
pbest=zeros(1,SearchAgents_no);
pbest_pos=zeros(dim,SearchAgents_no);
nfes=0;
j=0:(1/(max_nfes-1)):1; % Learning Probability Curve Pc
j=j*3;
Pc=(0.05+((0.45).*(exp(j)-exp(j(1)))./(exp(j(max_nfes))-exp(j(1))))); 
%   plot(Pc);
%    xlabel({'Number of evaluations'});
%    ylabel('probability curve')
for i=1:SearchAgents_no %Calculate all the objective values first
         nfes=nfes+1;
        Fitness(1,i)=fobj(X(:,i)');
       pbest(1,i)=Fitness(1,i);
       pbest_pos(:,i)=X(:,i);
end
        [fitness_sorted I]=sort(pbest);
        sorted_population=pbest_pos(:,I);
        Food_fitness=fitness_sorted(1);
         Food_pos=sorted_population(:,1);
         cg_curve(nfes)=Food_fitness;   
while nfes<max_nfes
      r=(ub-lb)/8;
       w=0.9-nfes*((0.9-0.4)/max_nfes);
      n=25;
    for i=1:n
        dz=randperm(25);
        r1=dz(1);    
        r2=dz(2);
        if r1==i
            r1=dz(6);
        elseif r2==i
            r2=dz(6);
        end  
       if rand>=Pc(nfes)
        for j=1:dim
            
          X(j,i)=X(j,i)+0.9*(X(j,r1)-X(j,r2))+0.9*(Food_pos(j)-X(j,i));
        end
       else
           for j=1:dim
       X(j,i)=X(j,i)+0.9*(sorted_population(j,r1)-sorted_population(j,r2))+0.9*(Food_pos(j)-X(j,i));
            end
             
        end
        Flag4ub=X(:,i)>ub';
        Flag4lb=X(:,i)<lb';
        X(:,i)=(X(:,i).*(~(Flag4ub+Flag4lb)))+rand*ub'.*Flag4ub+rand*lb'.*Flag4lb;
        Fitness(1,i)=fobj(X(:,i)');
        nfes=nfes+1;
         if Fitness(1,i)<pbest(1,i)
             pbest(1,i)=Fitness(1,i);
             pbest_pos(:,i)=X(:,i);
         end       
         cg_curve(nfes)=Food_fitness; 
          if mod(nfes,(500))==0
       display(['At NFES ', num2str(nfes), ' the best fitness is ', num2str(Food_fitness)]);
          end
    end
    for i=(n+1):SearchAgents_no
        index=0;
        neighbours_no=0;
        clear Neighbours_DeltaX
        clear Neighbours_X
        for j=(n+1):SearchAgents_no
            q=ceil(rand*3);
            Dist2Enemy=abs(X(q,i)-X(q,j));
            if (all(Dist2Enemy<=r(q)) && all(Dist2Enemy~=0))
                index=index+1;
                neighbours_no=neighbours_no+1;
                Neighbours_DeltaX(:,index)=DeltaX(:,j);
                Neighbours_X(:,index)=X(:,j);
            end
        end
        if neighbours_no>1
%          DeltaX(:,i)=(Food_pos-X(:,i));
        % Seperation
        % Eq. (3.1)
        S=zeros(dim,1);
        if neighbours_no>1
            for k=1:neighbours_no
                S=S+(Neighbours_X(:,k)-X(:,i));
            end
            S=-S;
        else
            S=zeros(dim,1);
        end
        
        % Alignment
        % Eq. (3.2)
        if neighbours_no>1
            A=(sum(Neighbours_DeltaX')')/neighbours_no;
        else
            A=DeltaX(:,i);
        end  
        
        if neighbours_no>1
            C_temp=(sum(Neighbours_X')')/neighbours_no;
        else
            C_temp=X(:,i);
        end 
        C=C_temp-X(:,i);
        
                for j=1:dim
                    DeltaX(j,i)=w*DeltaX(j,i)+rand*A(j,1)+rand*C(j,1)+rand*S(j,1);
                    if DeltaX(j,i)>Delta_max(j)
                        DeltaX(j,i)=Delta_max(j);
                    end
                    if DeltaX(j,i)<-Delta_max(j)
                        DeltaX(j,i)=-Delta_max(j);
                    end
                    X(j,i)=X(j,i)+DeltaX(j,i);
                end

        else   
                 X(:,i)=X(:,i)+Levy(dim)'.*X(:,i);
                 DeltaX(:,i)=0;
        end
        Flag4ub=X(:,i)>ub';
        Flag4lb=X(:,i)<lb';
        X(:,i)=(X(:,i).*(~(Flag4ub+Flag4lb)))+rand*ub'.*Flag4ub+rand*lb'.*Flag4lb;
        Fitness(1,i)=fobj(X(:,i)');
        nfes=nfes+1;
         if Fitness(1,i)<pbest(1,i)
             pbest(1,i)=Fitness(1,i);
             pbest_pos(:,i)=X(:,i);
         end
         cg_curve(nfes)=Food_fitness;   
    end
        [fitness_sorted I]=sort(pbest);
        sorted_population=pbest_pos(:,I);
        Food_fitness=fitness_sorted(1);
         Food_pos=sorted_population(:,1);
         if mod(nfes,(1000))==0
       display(['At NFES ', num2str(nfes), ' the best fitness is ', num2str(Food_fitness)]);
         end
end
Best_score=Food_fitness;
Best_pos=Food_pos;