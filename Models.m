clear
%to generate number of active coordinators (number of business opp.)
% S curve
x = 1:1:200;
fx = zeros(1,length(x));
a = 2; b=215;

for i = 1: length(x)
    if x(i) < a
        fx(i) = 0;
    else if x(i)>= a && x(i)<= (a+b)/2
            fx(i) = 2*((x(i)-a)/(b-a))^2;
        else if  x(i) >= (a+b)/2 && x(i)<=b
                fx(i) = 1-2*((x(i)-b)/(b-a))^2;
            else if x(i)>=b
                    fx(i) = 1;
                end
            end
        end
    end
end


%y = sigmf(x,[1 8]);
%plot(x,fx)


%Supplier
% long-term or short-term
% type of technology provided,
% degree of centrality
% number of successfull projects
% relation with other coordinators or partners

%input data
Data= struct; 

Data.number_ofShortTermSuppliers = 50;
Data.number_ofLongTermSuppliers = 50;
Data.number_ofSuppliers = Data.number_ofShortTermSuppliers + Data.number_ofLongTermSuppliers;

Data.number_ofReltCoordinators = 10;
Data.number_ofStrcCoordinators = 10;
Data.number_ofCoordinators = Data.number_ofReltCoordinators + Data.number_ofStrcCoordinators;

Data.number_ofTechnologies = 15;

Data.techMatrix_CR = zeros(Data.number_ofTechnologies, Data.number_ofCoordinators);  % 0 - 1 matrix (1 if coordinator has that technology)
Data.techMatrix_Spl = zeros(Data.number_ofTechnologies, Data.number_ofSuppliers);  % 0 - 1 matrix (1 if supplier has that technology)
% each supplier has only one type of technology

Data.Resource_CR = zeros(Data.number_ofTechnologies, Data.number_ofCoordinators);
Data.Resource_Spl = zeros(Data.number_ofTechnologies, Data.number_ofSuppliers);  %resources in terms of person month

Data.Cost_Spl = zeros(Data.number_ofTechnologies, Data.number_ofSuppliers); % person month cost
Data.Cost_CR = zeros(Data.number_ofTechnologies, Data.number_ofCoordinators);

Data.coordinator_resource = 3;  %initial assumption 
Data.cost_coeff_crd = 75;       %initial assumption
count_tech = zeros(1,Data.number_ofTechnologies); %counter for each technology type
for i = 1: Data.number_ofCoordinators
    cr_indx = 0;
    while cr_indx == 0
        tech_ind = randi([1 Data.number_ofTechnologies],1,1);
        if count_tech(1,tech_ind) < 2
            Data.techMatrix_CR(tech_ind,i) = 1;
            Data.Resource_CR(tech_ind,i) = Data.coordinator_resource;
            Data.Cost_CR(tech_ind,i) = Data.cost_coeff_crd;
            cr_indx = 1;
            count_tech(1,tech_ind) = count_tech(1,tech_ind)+1;
        end
    end
end

Data.supplier_resource = 4;  %initial assumption
Data.cost_coeff_sup = 50;    %initial assumption
count_tech = zeros(1,Data.number_ofTechnologies); %counter for each technology type
for i = 1:Data.number_ofSuppliers
    sup_indx = 0;
    while sup_indx == 0
        tech_ind = randi([1 Data.number_ofTechnologies],1,1);
        if count_tech(1,tech_ind) < 7
            Data.techMatrix_Spl(tech_ind,i) = 1;
            Data.Resource_Spl(tech_ind,i) = Data.supplier_resource;
            Data.Cost_Spl(tech_ind,i) = Data.cost_coeff_sup;
            sup_indx = 1;
            count_tech(1,tech_ind) = count_tech(1,tech_ind)+1;
        end
    end
end

Relation_Matrix = zeros(Data.number_ofCoordinators, Data.number_ofSuppliers);
Cumulative_ProfitMatrix = zeros(Data.number_ofCoordinators, Data.number_ofSuppliers);

sim_time = 200; %simulation time

max_active_crd = 15; %upper bound on active coordinators
%TU_weight = [10, 10, 8, 8, 8, 7, 7, 7, 6, 5, 5, 4, 4, 3, 2];  %remove this

%TU uncertaint - Generate input for sim_time
TU_uncrt = zeros(1,sim_time);
for t = 1: sim_time
    if rand < 0.5
        TU_uncrt(1,t) = 1;   %low TU
    else
        TU_uncrt(1,t) = 2;   %high TU
    end
end
bnd_high = [5, 7];  %uniform distribution bounds for business technology
bnd_low = [2, 4];

Data.req_resources = 4;  %required resources for each technology
Data.quality_indx = 1;   %for minimum technology

Data.Beta = 100; %long_term supplier scale parameter

%monitoring cost parameters
%monitoring_cost_coeff for each monitoring strategy
Data.m_s = [5, 10, 15, 20; 10, 20, 30, 40]; %Line1 cost for existing suppliers, Line2 for new suppliers

Data.Pd_s = [0.25, 0.50, 0.70, 0.90];  %Probability of detection

% Simulation
for t = 1: sim_time
    
    % Generate Business opportunity
    BussOpp = struct;  %Business opportunity is generated as a struct
    Act_Crd = struct;
    if t ==1
        %number of active coordinators
        Act_Crd.active_coord = 3;  % first we have 3 active coordinators
    else
        %number of active coordinators
        active_coord = max(max_active_crd*fx(t),active_coord);  
        Act_Crd.active_coord = round(active_coord);
    end
    
    %choose active coordinators
    Act_Crd.act_crd_indx = datasample(1:Data.number_ofCoordinators,Act_Crd.active_coord,'Replace',false); %randi([1 number_ofCoordinators],act_crd,1);
    
    tu_num = zeros(Act_Crd.active_coord,1);   %number of required technology for each business op
    %profit_margin = zeros(1,act_crd); %generate profit margin 0 < pm < 1 per person month
    req_technologies = [];
    for na = 1:Act_Crd.active_coord
        
        if TU_uncrt(1,t) == 1  %if TU is low
            tu_num(na) = randi(bnd_low,1,1);
        else
            tu_num(na) = randi(bnd_high,1,1);
        end
        BussOpp(na).Techindex = randsample(Data.number_ofTechnologies,tu_num(na));
        BussOpp(na).Techrequired = Data.req_resources*ones(tu_num(na),1);
        BussOpp(na).quality = Data.quality_indx*Data.req_resources*ones(tu_num(na),1);
        BussOpp(na).coordinator = Act_Crd.act_crd_indx(na);
        BussOpp(na).profit_margin = 0.2; %this can be randomly assigned
        
        %    %min resources for each technology  uniform between [10,40]
        %     min_res_tech = zeros(1,tu_num(na));
        %     for i = 1:tu_num(na)
        %         min_res_tech(i) = randi([1 5],1,1); %number of resources required for each technology, generated uniformly
        %     end
        req_technologies = [req_technologies BussOpp(na).Techindex'];
        
    end
    
    %%%%%%%%%%%%% bidding stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%Supplier Problem
    [Offered_coordinators, Offered_resources] = BiddingStagePartner(req_technologies,...
        Cumulative_ProfitMatrix, Data, BussOpp, Act_Crd);
    
    %%Coordinator Problem
    Selected_partners = PartnerSelection(BussOpp, Data, Act_Crd, Offered_coordinators,...
    Offered_resources, Cumulative_ProfitMatrix);

   %%% implementation stage
    
end
