

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end

Positions=Positions';
