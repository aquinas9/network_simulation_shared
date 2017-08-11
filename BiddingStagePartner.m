%%%% Bidding Stage
function [Offered_coordinators, Offered_resources] = BiddingStagePartner(req_technologies,...
        Cumulative_ProfitMatrix, Data, BussOpp, Act_Crd)

number_ofSuppliers = Data.number_ofSuppliers;
number_ofCoordinators = Data.number_ofCoordinators;
number_ofShortTermSuppliers = Data.number_ofShortTermSuppliers;
techMatrix_Spl = Data.techMatrix_Spl;
Resource_Spl = Data.Resource_Spl;
Cost_Spl = Data.Cost_Spl; 
Scale_param = Data.Beta;
num_ActCoordinator = Act_Crd.active_coord;
Act_Coordinator_indx = Act_Crd.act_crd_indx;

required_technologies = unique(req_technologies);
    Offered_coordinators = zeros(number_ofSuppliers, number_ofCoordinators);
    Offered_resources = zeros(number_ofSuppliers, 1);
    All_profits = sum(Cumulative_ProfitMatrix,2);
    Market_avg = sum(All_profits)/number_ofCoordinators;
    for sp = 1:number_ofSuppliers
        sup_tech = find(techMatrix_Spl(:,sp) == 1);  %supplied technology
        if any(sup_tech==required_technologies)
            %find business opportunity
            if sp <= number_ofShortTermSuppliers
                profit = zeros(1,num_ActCoordinator);
                for ac = 1:num_ActCoordinator
                    if any(sup_tech==BussOpp(ac).Techindex)
                        profit(ac) = BussOpp(ac).profit_margin*Cost_Spl(sup_tech,sp)*Resource_Spl(sup_tech,sp);
                    end
                end
                crd_indices = find(profit==max(profit));
                indx = randi(numel(crd_indices));
                offered_crd = crd_indices(indx);
                Offered_coordinators(sp,Act_Coordinator_indx(offered_crd))=1;
                Offered_resources(sp) = Resource_Spl(sup_tech, sp);
            else  %long-term supplier
                if Market_avg > 0
                    profit = zeros(1,num_ActCoordinator);
                    for ac = 1:num_ActCoordinator
                        if any(sup_tech==BussOpp(ac).Techindex)
                            profit(ac) = BussOpp(ac).profit_margin*Cost_Spl(sup_tech,sp)*Resource_Spl(sup_tech,sp)+...
                                Scale_param*(All_profits(Act_Coordinator_indx(ac))/Market_avg);
                        end
                    end
                else
                    profit = zeros(1,num_ActCoordinator);
                    for ac = 1:num_ActCoordinator
                        if any(sup_tech==BussOpp(ac).Techindex)
                            profit(ac) = BussOpp(ac).profit_margin*Cost_Spl(sup_tech,sp)*Resource_Spl(sup_tech,sp);
                        end
                    end
                end
                
                crd_indices = find(profit==max(profit));
                indx = randi(numel(crd_indices));
                offered_crd = crd_indices(indx);
                Offered_coordinators(sp,Act_Coordinator_indx(offered_crd))=1;
                Offered_resources(sp) = Resource_Spl(sup_tech, sp);
            end
        end
    end