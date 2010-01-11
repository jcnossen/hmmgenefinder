% HMM class for merging HMM objects.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko, Jelmer Cnossen, Orr Shomroni
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Merge < HMM
    methods
        % Class constructor.
        % Input:
        %       modelA - first model to merge
        %       modelB - second model to merge
        % Output:
        %       obj    - HMM object representing merged model.
        function [obj] = HMM_Merge(modelA, modelB)
            obj.StateCount = modelA.StateCount;
            obj.Trans = modelA.Trans;
            obj.Emit = modelA.Emit;
            obj.StateNames = modelA.StateNames;
            
            % Now to extend it with modelB
            bN = modelB.StateCount;
            ids = zeros(1, bN);
            for i = 1:bN
                AstateID = obj.get_id_by_state_name(modelB.StateNames{i});
                BstateID = modelB.get_id_by_state_name(modelB.StateNames{i});
                if (AstateID == 0)
                    % no such node in A - just add one
                    AstateID = obj.add_state(modelB.StateNames{i}, modelB.Emit(BstateID, :));
                else
                    fprintf('[!] Duplicate state found: %s (%i)\n', modelB.StateNames{i}, AstateID);
                end
                ids(BstateID) = AstateID;
            end
            
            % Now add all edges
            for i = 1:bN
                for j = 1:bN
                    p = modelB.Trans(i, j);
                    if (p > eps)
                        if (obj.Trans(ids(i), ids(j)) > 0 && obj.Trans(ids(i), ids(j)) ~= p)
                            fprintf('[!] Duplicate edge found: %s (%i) -> %s (%i), overwriting probability %.5f -> %.5f\n', obj.StateNames{ids(i)}, ids(i), obj.StateNames{ids(j)}, ids(j), obj.Trans(ids(i), ids(j)), p );
                        end
                        obj.add_edge(ids(i), ids(j), p);
                    end
                end
            end
        end
    end
end