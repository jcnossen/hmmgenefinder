% Builds an HMM for our DBDM assignment.
% Input arguments:
%      <seq>       = training set for the algorithm, as output by
%                    split_ecoli command.
%      <modelType> = type of model to use, possible types are 'dumb',
%                    'full'
%      [*prob*]    = insertion and deletion probabilities used by Genic
%                    model.
% Output arguments:
%      [resHMM]    = object of class HMM, which contains all info required
%                    to define an HMM.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [resHMM] = build_hmm(seq, modelType, insertion_probability, deletion_probability) 
    % Making life easy
    if (nargin < 2)
        error('Not enough parameters, consult help.');
    end
    
    modelType = find(strcmpi(modelType, {'dumb', 'full'}), 1);
    if (isempty(modelType))
        error('Invalid model type specified');
    end
    
    if (nargin < 4)
        insertion_probability = 10 ^ -8;
        deletion_probability = 10 ^ -8;
        fprintf('[!] No insertion or deletion probabilities defined.\n');
        fprintf('[i] Assuming Pins = %e, Pdel = %e\n', insertion_probability, deletion_probability);
    end

    resHMM = [];
    
    genic_part = HMM_Genic(seq, 'complex', insertion_probability, deletion_probability);
    if (modelType == 2) % full
        intergenic_short_part = HMM_Intergenic_Short(seq, true);
        intergenic_long_part = HMM_Intergenic_Long(seq, true);
        
        % Merge all into 1 model
        resHMM = HMM_Merge(HMM_Merge(genic_part, intergenic_short_part), intergenic_long_part);
        
        % Link intergenic -> start codons
        resHMM.add_edge('intergenic_short_end', 'start_codons_AGT', 1);
        resHMM.add_edge('intergenic_long_end', 'start_codons_AGT', 1);
        
        % Link stop codons -> intergenic
        intID = resHMM.add_state('intergenic_begin', [0 0 0 0]);
        
        pShort = intergenic_short_part.GeneCount / (intergenic_short_part.GeneCount + intergenic_long_part.GeneCount);
        pLong = intergenic_long_part.GeneCount / (intergenic_short_part.GeneCount + intergenic_long_part.GeneCount);
        resHMM.add_edge(intID, 'intergenic_short_begin', pShort);
        resHMM.add_edge(intID, 'intergenic_long_begin', pLong);
        resHMM.add_edge('stop_TAA_TGA_3', intID, 1);
        resHMM.add_edge('stop_TAG_3', intID, 1);
    end
    
    if (modelType == 1) % dumb
        %intergenic_part = HMM_Intergenic_Dumb(seq, false);
        intergenic_part = HMM_Intergenic_Dumb(seq, true);
        
        % Merge all into 1 model
        resHMM = HMM_Merge(genic_part, intergenic_part, [], [true true]);
    end
end