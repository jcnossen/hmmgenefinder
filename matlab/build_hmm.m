% Builds an HMM for our DBDM assignment.
% Input arguments:
%      <seq>        = training set for the algorithm, as output by
%                     split_ecoli command.
%      <modelType>  = type of model to use, possible types are 'dumb',
%                     'full'
%      [ins_p]      = probability of inserting a nucleotide in genic model
%      [del_p]       = proability of deleting a nucleotide in genic model
%      [trainGenic] = boolean switch determining whether genic model must
%                     be trained. Default - false.
%      [trainInter] = boolean switch determining whether intergenic model
%                     should be trained. Default - false.
% Output arguments:
%      [resHMM]    = object of class HMM, which contains all info required
%                    to define an HMM.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [resHMM] = build_hmm(seq, modelType, ins_p, del_p, trainGenic, trainInter) 
    % Making life easy
    if (nargin < 2)
        error('Not enough parameters, consult help.');
    end
    
    modelType = find(strcmpi(modelType, {'dumb', 'full'}), 1);
    if (isempty(modelType))
        error('Invalid model type specified');
    end
    
    if (nargin < 4)
        ins_p = 10 ^ -8;
        del_p = 10 ^ -8;
        fprintf('[!] No insertion or deletion probabilities defined.\n');
        fprintf('[i] Assuming Pins = %e, Pdel = %e\n', ins_p, del_p);
    end
    
    if (nargin < 6)
        fprintf('[!] Warning: model training is disabled by default.\n');
        trainGenic = false;
        trainInter = false;
    end

    resHMM = [];
    
    fprintf('[i] Build: creating genic model.\n');
    genic = HMM_Genic(seq, 'complex', ins_p, del_p, trainGenic);
    
    if (modelType == 2) % full
        fprintf('[i] Build: creating short intergenic model.\n');
        intergenic_short = HMM_Intergenic_Short(seq, trainInter);
        fprintf('[i] Build: creating long intergenic model.\n');
        intergenic_long = HMM_Intergenic_Long(seq, trainInter);
        
        % Probabilities
        shortCount = intergenic_short.GeneCount;
        longCount = intergenic_long.GeneCount;
        both = shortCount + longCount;
        
        % Short
        IDs = [intergenic_short.get_id_by_state_name('stop_TAG_3') intergenic_short.get_id_by_state_name('stop_TAA_TGA_3')];
        intergenic_short.Trans(IDs, :) = intergenic_short.Trans(IDs, :) .* (shortCount / both);
        % Long
        IDs = [intergenic_long.get_id_by_state_name('stop_TAG_3') intergenic_long.get_id_by_state_name('stop_TAA_TGA_3')];
        intergenic_long.Trans(IDs, :) = intergenic_long.Trans(IDs, :) .* (longCount / both);
        
        
        % Create intergenic model
        fprintf('[i] Build: merging long and short intergenic models.\n');
        intergenic = HMM_Merge(intergenic_long, intergenic_short, [], [false false]);
        fprintf('[i] Build: merging genic and intergenic models.\n');
        resHMM = HMM_Merge(genic, intergenic);
    end
    
    if (modelType == 1) % dumb
        fprintf('[i] Build: creating intergenic model.\n');
        intergenic = HMM_Intergenic_Dumb(seq, trainInter);
        
        % Merge all into 1 model
        fprintf('[i] Build: merging genic and intergenic models.\n');
        resHMM = HMM_Merge(genic, intergenic, [], [true true]);
    end
end