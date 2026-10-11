% Loading data from excel
data = readtable("Data for G-DEA three groups.xlsx");

% Get the DMU name column and separate it
DMU = data.DMU;

% Creating BOLECREVA
BOLECREVA = table2array(data(:, ...
    ["Inv_IE","Inv_OE","R","TA","ES","SS","GS"]));

save BOLECREVA BOLECREVA DMU                                                                                                            