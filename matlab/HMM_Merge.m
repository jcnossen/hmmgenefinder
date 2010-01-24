% HMM class for merging HMM objects.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Merge < HMM
    methods
        % Class constructor.
        % Input:
        %       <modelA>              - first model to merge
        %       <modelB>              - second model to merge
        %       [stateConflicts]      - array with boolean values defining
        %                               which of two concurent states to
        %                               keep (true for modelA, false for
        %                               modelB)
        %       [transitionConflicts] - array with boolean values defining
        %                               which of two concurent transitions
        %                               to keep (true for modelA, false for
        %                               modelB). This is used only when
        %                               both states have non-zero
        %                               transition probability
        % Output:
        %       [obj]    - HMM object representing merged model.
        function [obj] = HMM_Merge(modelA, modelB, stateConflicts, transitionConflicts)
            if (nargin < 2)
                error('Not enough input arguments. Consult help');
            end
            if (nargin < 3)
                stateConflicts = [];
            end
            if (nargin < 4)
                transitionConflicts = [];
            end
            
            obj.StateCount = modelA.StateCount;
            obj.Trans = modelA.Trans;
            obj.Emit = modelA.Emit;
            obj.StateNames = modelA.StateNames;
            
            % Now to extend it with modelB
            stateDecision = 0;
            nStateConflicts = length(stateConflicts);
            bN = modelB.StateCount;
            ids = zeros(1, bN);
            for i = 1:bN
                AstateID = obj.get_id_by_state_name(modelB.StateNames{i});
                BstateID = modelB.get_id_by_state_name(modelB.StateNames{i});
                if (AstateID == 0)
                    % no such node in A - just add one
                    AstateID = obj.add_state(modelB.StateNames{i}, modelB.Emit(BstateID, :));
                else
                    stateDecision = stateDecision + 1;
                    if (stateDecision > nStateConflicts || stateConflicts(stateDecision))
                        fprintf('[!] Duplicate state found: %s (%i) -> keeping emission probabilities.\n', modelB.StateNames{i}, AstateID);
                    else%if (stateDecision <= nStateConflicts && ~stateConflicts(stateDecision))
                        obj.Emit(AstateID, :) = modelB.Emit(BstateID, :);
                        fprintf('[!] Duplicate state found: %s (%i) -> overwriting emission probabilities.\n', modelB.StateNames{i}, AstateID);
                    end
                end
                ids(BstateID) = AstateID;
            end
            
            transitionDecision = 0;
            nTransitionConflicts = length(transitionConflicts);
            % Now add all edges
            for i = 1:bN
                for j = 1:bN
                    p = modelB.Trans(i, j);
                    if (p > eps)
                        addEdge = true;
                        if (obj.Trans(ids(i), ids(j)) > 0 && obj.Trans(ids(i), ids(j)) ~= p)
                            transitionDecision = transitionDecision + 1;
                            if (transitionDecision > nTransitionConflicts || transitionConflicts(transitionDecision))
                                fprintf('[!] Duplicate edge found: %s (%i) -> %s (%i), keeping probability %.5f\n', obj.StateNames{ids(i)}, ids(i), obj.StateNames{ids(j)}, ids(j), obj.Trans(ids(i), ids(j)));
                                addEdge = false;
                            else
                                fprintf('[!] Duplicate edge found: %s (%i) -> %s (%i), overwriting probability %.5f -> %.5f\n', obj.StateNames{ids(i)}, ids(i), obj.StateNames{ids(j)}, ids(j), obj.Trans(ids(i), ids(j)), p );
                            end
                        end
                        if (addEdge)
                            obj.add_edge(ids(i), ids(j), p);
                        end
                    end
                end
            end
        end
    end
end