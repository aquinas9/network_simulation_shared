   
function z_j  = PartnerSelection(BussOpp, Data, Act_Crd, Offered_coordinators,...
    Offered_resources, Cumulative_ProfitMatrix)

number_ofReltCoordinators = Data.number_ofReltCoordinators;
num_ActCoordinator = Act_Crd.active_coord;
Act_Coordinator_indx = Act_Crd.act_crd_indx;
m_s = Data.m_s;
Cost_Spl = Data.Cost_Spl; 
Scale_param = Data.Beta;
Resource_CR = Data.Resource_CR;
Cost_CR = Data.Cost_CR;


%%Cooridnators Problem %%%%%%%%%%%%%%
    for na = 1:num_ActCoordinator
        
        if  Act_Coordinator_indx(na) <= number_ofReltCoordinators  %coordinator is a relational coordinator
            list_ofPossibleSuppliers = find(Offered_coordinators(:, Act_Coordinator_indx(na)) > 0);
            Resources_forCrd = Offered_resources(list_ofPossibleSuppliers);
            Delta = zeros(1, length(list_ofPossibleSuppliers));
            Qual = zeros(1, length(list_ofPossibleSuppliers));
              
            %f(x,h) = Delta*x_ij measure of the quality of the supplier?s performance
            cum_profit = Cumulative_ProfitMatrix(na,list_ofPossibleSuppliers);
            avg_profit = mean(cum_profit);
            if avg_profit == 0  % there is no previous experience
               Delta(1,:) = 1;
               Qual(1,:)=Delta(1,:).*Resources_forCrd';
            else
               for i=1:length(list_ofPossibleSuppliers)
                   if cum_profit(list_ofPossibleSuppliers(i)) >= avg_profit
                      Delta(1,i) = 0.75;
                   else
                      Delta(1,i) = 0.50;
                   end
               end
              indx_max_profit = find(cum_profit(1,:) == max(cum_profit(1,:)));
              Delta(na,indx_max_profit) = 1;
              Qual(1,:)=Delta(1,:).*Resources_forCrd';
            end
           
            %Monitoring strategy and monitoring cost
            %each strategy is a decision 
            MC = zeros(length(list_ofPossibleSuppliers), length(m_s(1,:))); %Monitoring cost matrix columns are strategies       
            for j=1:length(list_ofPossibleSuppliers)
                if Cumulative_ProfitMatrix(Act_Coordinator_indx(na),list_ofPossibleSuppliers(j)) > 0 %coordinator knows the supplier  
                    MC(j,:) = m_s(1,:)*Resources_forCrd(j);   %monitoring cost for each strategy
                else  %coordinator does not know the supplier
                    MC(j,:) = m_s(2,:)*Resources_forCrd(j);   %monitoring cost for each strategy
                end
            end
            
            %optimisation model of coordinator na
            %decision variables which partners to choose, 
            z_j = zeros(length(list_ofPossibleSuppliers),1);
            %which monitorin strategy to choose and 
            y_js = zeros(length(list_ofPossibleSuppliers),length(m_s));
            %how many resourses to allocate to technology i
            u_i = zeros(1,1);
            total_dv = length(list_ofPossibleSuppliers) + length(list_ofPossibleSuppliers)*length(m_s) + 1;
            
            %Construct objective function
            Cost_ui =  sum(Cost_CR(:,Act_Coordinator_indx(na))); %cost of technology i for coordinator
            Cost_coeff_zj = zeros(length(list_ofPossibleSuppliers),1);
            Preference_zj = zeros(length(list_ofPossibleSuppliers),1);
            Cost_coeff_yjs = MC;
            for i=1:length(list_ofPossibleSuppliers)
                Cost_coeff_zj(i) = BussOpp(na).profit_margin*sum(Cost_Spl(:,list_ofPossibleSuppliers(i)));
                Preference_zj(i) = Scale_param*Delta(1,i);
            end
            Coeff_zj = Cost_coeff_zj - Preference_zj;
            
           %Construct Constraints
           %Constraint 1 ==>   u_i <= b_i   %coordinator only supplies one technology
           C_1 = zeros(1, total_dv); C_1(1,1) = 1;  Cb_1 = sum(Resource_CR(:,Act_Coordinator_indx(na)));
           
           %Constraint 2 ==>   
           % for each strategy solve the implementation problem and obtain how much partner will defect
           
           %%% implementation stage for each strategy
           for st = 1:length(m_s)
               
           end
           
           
            
        else   %coordinator is a structural coordinator
            
            
        end
    end