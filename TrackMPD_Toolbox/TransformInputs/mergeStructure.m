function [resultStruct] = mergeStructure(mainStruct,struct2merge)
fields2merge = fields(struct2merge);
for ifieldsIn = fields2merge'
    moreLevels = isstruct(struct2merge.(ifieldsIn{1}));
    if moreLevels 
        [valueInStructLeveli] = mergeStructure(mainStruct.(ifieldsIn{1}),struct2merge.(ifieldsIn{1}));
        mainStruct.(ifieldsIn{1}) = valueInStructLeveli;
    else
        mainStruct.(ifieldsIn{1}) = struct2merge.(ifieldsIn{1});
    end
end
resultStruct = mainStruct;
end