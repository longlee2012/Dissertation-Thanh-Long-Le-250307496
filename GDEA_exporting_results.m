Efficiency = table(dea_ahp_results);
writetable(Efficiency, 'GDEA_efficiency_3_groups_with_restrictions.xlsx');
Weights = array2table(weights_for_each_alternative', ...
    'VariableNames', ...
    {'Inv_IE','Inv_OE','R','TA','ES','SS','GS'});
Weights = addvars(Weights, DMU, 'Before', 1);
writetable(Weights, 'GDEA_weights_3_groups_with_restrictions.xlsx');
