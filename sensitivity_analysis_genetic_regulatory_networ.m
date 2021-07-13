%this script is to perform sensitivity analysis on genetic transcriptions,
%with the goal of identifying the influential genetic regulators in a
%genetic regulatory network, hoping to decrease the model order. 

%in this example, the experiment is run 6 times, with the same genetic
%candidates in the analysis, for robustness and benchmarking. 

%the following script analyses the sensitivity by regression analysis
%for ease of understanding and reproducibility, the code is basically
%repeated six times. In practice, you can of course put the reptitions in
%loops. 

x = readtable('out1.csv')
input_gene = x(:,(1:14))
output_gene = x(:,15)
a = table2array(input_gene(:,1))
b= table2array(input_gene(:,2))
c= table2array(input_gene(:,3))
d= table2array(input_gene(:,4))
e= table2array(input_gene(:,5))
f= table2array(input_gene(:,6))
g= table2array(input_gene(:,7))
h= table2array(input_gene(:,8))
i= table2array(input_gene(:,9))
j= table2array(input_gene(:,10))
l= table2array(input_gene(:,11))
m= table2array(input_gene(:,12))
n= table2array(input_gene(:,13))
o = table2array(input_gene(:,14))

%defining input benchmark
input_gene = [a b c d e f g h i j l m n o];
input_gene = normalize(input_gene);
%defining output benchmark
%the input and output definitions could change, depending on how one would
%define regulators 
output_gene = table2array(output_gene);
runorder = randperm(10);
Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
    input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','LBD37','HHO3','HHO2','ZFP4','BZIP53','NAP57','AT4G12750','BZO2H1','XERICO','AT2G43500','MYB42','MYB61','VND1','AT1G03040','AT1G08090_1'});
% Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
%     input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'});

modelspec = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';

%fitting linear model, with the goal of regression
mdl = fitlm(Expmt,modelspec);
aaa = find(mdl.Coefficients.Estimate)
name111 = mdl.CoefficientNames(aaa)
model10 = mdl.Coefficients.Estimate(aaa);
figure()
h = bar(mdl.Coefficients.Estimate(aaa));
set(h,'facecolor',[0.8 0.8 0.9])
legend('Coefficient')
set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
set(gca,'xticklabel',[name111(1) name111(2) name111(3) name111(4) name111(5) name111(6) name111(7) name111(8) name111(9) name111(10)])
ylabel('Input Genes')
xlabel('Normalized Coefficient')
title('Quadratic Model Coefficients')

Experiment_11 = table(name111',model10);
Experiment_11.Properties.VariableNames = {'Gene' 'Coefficient'};
writetable(Experiment_11,'Experiment 11.xlsx')

% 
% plotSlice(mdl)


%%
x2 = readtable('out2.csv')

input_gene2 = x2(:,(1:14))
output_gene2 = x2(:,15)
a2 = table2array(input_gene2(:,1))
b2= table2array(input_gene2(:,2))
c2= table2array(input_gene2(:,3))
d2= table2array(input_gene2(:,4))
e2= table2array(input_gene2(:,5))
f2= table2array(input_gene2(:,6))
g2= table2array(input_gene2(:,7))
h2= table2array(input_gene2(:,8))
i2= table2array(input_gene2(:,9))
j2= table2array(input_gene2(:,10))
l2= table2array(input_gene2(:,11))
m2= table2array(input_gene2(:,12))
n2= table2array(input_gene2(:,13))
o2 = table2array(input_gene2(:,14))

input_gene2 = [a2 b2 c2 d2 e2 f2 g2 h2 i2 j2 l2 m2 n2 o2];
input_gene2 = normalize(input_gene2);
output_gene2 = table2array(output_gene2);
% input_gene = table2array(input_gene);
% output_gene = table2array(output_gene);
runorder = randperm(10);
Expmt2 = table(runorder', input_gene2(:,1), input_gene2(:,2), input_gene2(:,3), input_gene2(:,4), input_gene2(:,5), input_gene2(:,6),input_gene2(:,7), input_gene2(:,8),...
    input_gene2(:,9),input_gene2(:,10),input_gene2(:,11),input_gene2(:,12),input_gene2(:,13),input_gene2(:,14), output_gene2(:,1),'VariableNames',{'RunNumber','LBD37','HHO3','HHO2','ZFP4','BZIP53','NAP57','AT4G12750','BZO2H1','XERICO','AT2G43500','MYB42','MYB61','VND1','AT1G03040','AT1G08090_1'});
% Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
%     input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'});
% mdl = fitlm(Expmt,'o~a*b*c*d*e*f*g*h*i*j*k*l*m*n');
% modelspec = 'AT1G08090_1~ AT5G67420 + AT5G67420^2 + AT1G25550 + AT1G25550^2 +AT1G68670+ AT1G68670^2 +AT1G66140+ AT1G66140^2 +AT3G62420+AT3G62420^2 +AT3G57150+AT3G57150^2 + AT4G12750+ AT4G12750^2+AT4G02640^2+AT4G02640 +AT2G04240+AT2G04240^2+AT2G43500+AT2G43500^2+AT4G12350+ AT4G12350^2+AT1G09540+AT1G09540^2+AT2G18060+AT2G18060^2 +AT1G03040+ AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040+AT5G67420^2+AT1G25550^2+AT1G68670^2+AT1G66140^2 +AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec2 = 'AT1G08090_1~LBD37*HHO3*HHO2*ZFP4*BZIP53*NAP57*AT4G12750*BZO2H1*XERICO*AT2G43500*MYB42*MYB61*VND1*AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';
modelspec2 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';
% modelspec2 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040';
% modelspec = 'AT1G08090_1~ AT5G67420 * AT5G67420^2 * AT1G25550 * AT1G25550^2 *AT1G68670*AT1G68670^2 *AT1G66140*AT1G66140^2 *AT3G62420*AT3G62420^2 *AT3G57150*AT3G57150^2 *AT4G12750*AT4G12750^2*AT4G02640^2*AT4G02640*AT2G04240*AT2G04240^2*AT2G43500*AT2G43500^2*AT4G12350*AT4G12350^2*AT1G09540*AT1G09540^2*AT2G18060*AT2G18060^2 *AT1G03040*AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550^2+AT1G68670^2+AT1G66140^2+AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550+AT1G68670+AT1G66140+AT3G62420+AT3G57150+AT4G12750+AT4G02640+AT2G04240+AT2G43500+AT4G12350+AT1G09540+AT2G18060+AT1G03040';

% mdl = fitlm(Expmt,'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040');
% mdl = fitlm(Expmt,modelspec);
mdl2 = fitlm(Expmt2,modelspec2);
bbb = find(mdl2.Coefficients.Estimate)
name222 = mdl2.CoefficientNames(bbb)
model20 = mdl2.Coefficients.Estimate(bbb);
figure()
h2 = bar(mdl2.Coefficients.Estimate(bbb));
set(h2,'facecolor',[0.8 0.8 0.9])
legend('Coefficient')
set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
set(gca,'xticklabel',[name222(1) name222(2) name222(3) name222(4) name222(5) name222(6) name222(7) name222(8) name222(9) name222(10)])
ylabel('Input Genes')
xlabel('Normalized Coefficient')
title('Quadratic Model Coefficients')

Experiment_22 = table(name222',model20);
Experiment_22.Properties.VariableNames = {'Gene' 'Coefficient'};
writetable(Experiment_22,'Experiment 22.xlsx')

% plotSlice(mdl2)

%%
x3 = readtable('out3.csv')

input_gene3 = x3(:,(1:14))
output_gene3 = x3(:,15)
a3 = table2array(input_gene3(:,1))
b3= table2array(input_gene3(:,2))
c3= table2array(input_gene3(:,3))
d3= table2array(input_gene3(:,4))
e3= table2array(input_gene3(:,5))
f3= table2array(input_gene3(:,6))
g3= table2array(input_gene3(:,7))
h3= table2array(input_gene3(:,8))
i3= table2array(input_gene3(:,9))
j3= table2array(input_gene3(:,10))
l3= table2array(input_gene3(:,11))
m3= table2array(input_gene3(:,12))
n3= table2array(input_gene3(:,13))
o3 = table2array(input_gene3(:,14))

input_gene3 = [a3 b3 c3 d3 e3 f3 g3 h3 i3 j3 l3 m3 n3 o3];
input_gene3 = normalize(input_gene3);
output_gene3 = table2array(output_gene3);
% input_gene = table2array(input_gene);
% output_gene = table2array(output_gene);
runorder = randperm(10);
Expmt3 = table(runorder', input_gene3(:,1), input_gene3(:,2), input_gene3(:,3), input_gene3(:,4), input_gene3(:,5), input_gene3(:,6),input_gene3(:,7), input_gene3(:,8),...
    input_gene3(:,9),input_gene3(:,10),input_gene3(:,11),input_gene3(:,12),input_gene3(:,13),input_gene3(:,14), output_gene3(:,1),'VariableNames',{'RunNumber','LBD37','HHO3','HHO2','ZFP4','BZIP53','NAP57','AT4G12750','BZO2H1','XERICO','AT2G43500','MYB42','MYB61','VND1','AT1G03040','AT1G08090_1'});
% Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
%     input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'});
% mdl = fitlm(Expmt,'o~a*b*c*d*e*f*g*h*i*j*k*l*m*n');
% modelspec = 'AT1G08090_1~ AT5G67420 + AT5G67420^2 + AT1G25550 + AT1G25550^2 +AT1G68670+ AT1G68670^2 +AT1G66140+ AT1G66140^2 +AT3G62420+AT3G62420^2 +AT3G57150+AT3G57150^2 + AT4G12750+ AT4G12750^2+AT4G02640^2+AT4G02640 +AT2G04240+AT2G04240^2+AT2G43500+AT2G43500^2+AT4G12350+ AT4G12350^2+AT1G09540+AT1G09540^2+AT2G18060+AT2G18060^2 +AT1G03040+ AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040+AT5G67420^2+AT1G25550^2+AT1G68670^2+AT1G66140^2 +AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec3 = 'AT1G08090_1~LBD37*HHO3*HHO2*ZFP4*BZIP53*NAP57*AT4G12750*BZO2H1*XERICO*AT2G43500*MYB42*MYB61*VND1*AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';
% modelspec3 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040';
modelspec3 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';

% modelspec = 'AT1G08090_1~ AT5G67420 * AT5G67420^2 * AT1G25550 * AT1G25550^2 *AT1G68670*AT1G68670^2 *AT1G66140*AT1G66140^2 *AT3G62420*AT3G62420^2 *AT3G57150*AT3G57150^2 *AT4G12750*AT4G12750^2*AT4G02640^2*AT4G02640*AT2G04240*AT2G04240^2*AT2G43500*AT2G43500^2*AT4G12350*AT4G12350^2*AT1G09540*AT1G09540^2*AT2G18060*AT2G18060^2 *AT1G03040*AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550^2+AT1G68670^2+AT1G66140^2+AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550+AT1G68670+AT1G66140+AT3G62420+AT3G57150+AT4G12750+AT4G02640+AT2G04240+AT2G43500+AT4G12350+AT1G09540+AT2G18060+AT1G03040';

% mdl = fitlm(Expmt,'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040');
% mdl = fitlm(Expmt,modelspec);
mdl3 = fitlm(Expmt3,modelspec3);
ccc = find(mdl3.Coefficients.Estimate)
name333 = mdl3.CoefficientNames(ccc)
model30 = mdl3.Coefficients.Estimate(ccc);
figure()
h3 = bar(mdl3.Coefficients.Estimate(ccc));
set(h3,'facecolor',[0.8 0.8 0.9])
legend('Coefficient')
set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
set(gca,'xticklabel',[name333(1) name333(2) name333(3) name333(4) name333(5) name333(6) name333(7) name333(8) name333(9) name333(10)])
ylabel('Input Genes')
xlabel('Normalized Coefficient')
title('Quadratic Model Coefficients')

Experiment_33 = table(name333',model30);
Experiment_33.Properties.VariableNames = {'Gene' 'Coefficient'};
writetable(Experiment_33,'Experiment 33.xlsx')
% 
% plotSlice(mdl3)

%%
x4 = readtable('out4.csv')

input_gene4 = x4(:,(1:14))
output_gene4 = x4(:,15)
a4 = table2array(input_gene4(:,1))
b4= table2array(input_gene4(:,2))
c4= table2array(input_gene4(:,3))
d4= table2array(input_gene4(:,4))
e4= table2array(input_gene4(:,5))
f4= table2array(input_gene4(:,6))
g4= table2array(input_gene4(:,7))
h4= table2array(input_gene4(:,8))
i4= table2array(input_gene4(:,9))
j4= table2array(input_gene4(:,10))
l4= table2array(input_gene4(:,11))
m4= table2array(input_gene4(:,12))
n4= table2array(input_gene4(:,13))
o4 = table2array(input_gene4(:,14))

input_gene4 = [a4 b4 c4 d4 e4 f4 g4 h4 i4 j4 l4 m4 n4 o4];
input_gene4 = normalize(input_gene4);
output_gene4 = table2array(output_gene4);
% input_gene = table2array(input_gene);
% output_gene = table2array(output_gene);
runorder = randperm(10);
Expmt4 = table(runorder', input_gene4(:,1), input_gene4(:,2), input_gene4(:,3), input_gene4(:,4), input_gene4(:,5), input_gene4(:,6),input_gene4(:,7), input_gene4(:,8),...
    input_gene4(:,9),input_gene4(:,10),input_gene4(:,11),input_gene4(:,12),input_gene4(:,13),input_gene4(:,14), output_gene4(:,1),'VariableNames',{'RunNumber','LBD37','HHO3','HHO2','ZFP4','BZIP53','NAP57','AT4G12750','BZO2H1','XERICO','AT2G43500','MYB42','MYB61','VND1','AT1G03040','AT1G08090_1'});
% Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
%     input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'});
% mdl = fitlm(Expmt,'o~a*b*c*d*e*f*g*h*i*j*k*l*m*n');
% modelspec = 'AT1G08090_1~ AT5G67420 + AT5G67420^2 + AT1G25550 + AT1G25550^2 +AT1G68670+ AT1G68670^2 +AT1G66140+ AT1G66140^2 +AT3G62420+AT3G62420^2 +AT3G57150+AT3G57150^2 + AT4G12750+ AT4G12750^2+AT4G02640^2+AT4G02640 +AT2G04240+AT2G04240^2+AT2G43500+AT2G43500^2+AT4G12350+ AT4G12350^2+AT1G09540+AT1G09540^2+AT2G18060+AT2G18060^2 +AT1G03040+ AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040+AT5G67420^2+AT1G25550^2+AT1G68670^2+AT1G66140^2 +AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec4 = 'AT1G08090_1~LBD37*HHO3*HHO2*ZFP4*BZIP53*NAP57*AT4G12750*BZO2H1*XERICO*AT2G43500*MYB42*MYB61*VND1*AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';
% modelspec4 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040';
modelspec4 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';

% modelspec = 'AT1G08090_1~ AT5G67420 * AT5G67420^2 * AT1G25550 * AT1G25550^2 *AT1G68670*AT1G68670^2 *AT1G66140*AT1G66140^2 *AT3G62420*AT3G62420^2 *AT3G57150*AT3G57150^2 *AT4G12750*AT4G12750^2*AT4G02640^2*AT4G02640*AT2G04240*AT2G04240^2*AT2G43500*AT2G43500^2*AT4G12350*AT4G12350^2*AT1G09540*AT1G09540^2*AT2G18060*AT2G18060^2 *AT1G03040*AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550^2+AT1G68670^2+AT1G66140^2+AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550+AT1G68670+AT1G66140+AT3G62420+AT3G57150+AT4G12750+AT4G02640+AT2G04240+AT2G43500+AT4G12350+AT1G09540+AT2G18060+AT1G03040';

% mdl = fitlm(Expmt,'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040');
% mdl = fitlm(Expmt,modelspec);
mdl4 = fitlm(Expmt4,modelspec4);
ddd = find(mdl4.Coefficients.Estimate)
name444 = mdl4.CoefficientNames(ddd)
model40 = mdl4.Coefficients.Estimate(ddd);
figure()
h4 = bar(mdl4.Coefficients.Estimate(ddd));
set(h4,'facecolor',[0.8 0.8 0.9])
legend('Coefficient')
set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
set(gca,'xticklabel',[name444(1) name444(2) name444(3) name444(4) name444(5) name444(6) name444(7) name444(8) name444(9) name444(10)])
ylabel('Input Genes')
xlabel('Normalized Coefficient')
title('Quadratic Model Coefficients')

Experiment_44 = table(name444',model40);
Experiment_44.Properties.VariableNames = {'Gene' 'Coefficient'};
writetable(Experiment_44,'Experiment 44.xlsx')
% plotSlice(mdl4)

%%
x5 = readtable('out5.csv')

input_gene5 = x5(:,(1:14))
output_gene5 = x5(:,15)
a5 = table2array(input_gene5(:,1))
b5= table2array(input_gene5(:,2))
c5= table2array(input_gene5(:,3))
d5= table2array(input_gene5(:,4))
e5= table2array(input_gene5(:,5))
f5= table2array(input_gene5(:,6))
g5= table2array(input_gene5(:,7))
h5= table2array(input_gene5(:,8))
i5= table2array(input_gene5(:,9))
j5= table2array(input_gene5(:,10))
l5= table2array(input_gene5(:,11))
m5= table2array(input_gene5(:,12))
n5= table2array(input_gene5(:,13))
o5 = table2array(input_gene5(:,14))

input_gene5 = [a5 b5 c5 d5 e5 f5 g5 h5 i5 j5 l5 m5 n5 o5];
input_gene5 = normalize(input_gene5);
output_gene5 = table2array(output_gene5);
% input_gene = table2array(input_gene);
% output_gene = table2array(output_gene);
runorder = randperm(10);
Expmt5 = table(runorder', input_gene5(:,1), input_gene5(:,2), input_gene5(:,3), input_gene5(:,4), input_gene5(:,5), input_gene5(:,6),input_gene5(:,7), input_gene5(:,8),...
    input_gene5(:,9),input_gene5(:,10),input_gene5(:,11),input_gene5(:,12),input_gene5(:,13),input_gene5(:,14), output_gene5(:,1),'VariableNames',{'RunNumber','LBD37','HHO3','HHO2','ZFP4','BZIP53','NAP57','AT4G12750','BZO2H1','XERICO','AT2G43500','MYB42','MYB61','VND1','AT1G03040','AT1G08090_1'});
% Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
%     input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'});
% mdl = fitlm(Expmt,'o~a*b*c*d*e*f*g*h*i*j*k*l*m*n');
% modelspec = 'AT1G08090_1~ AT5G67420 + AT5G67420^2 + AT1G25550 + AT1G25550^2 +AT1G68670+ AT1G68670^2 +AT1G66140+ AT1G66140^2 +AT3G62420+AT3G62420^2 +AT3G57150+AT3G57150^2 + AT4G12750+ AT4G12750^2+AT4G02640^2+AT4G02640 +AT2G04240+AT2G04240^2+AT2G43500+AT2G43500^2+AT4G12350+ AT4G12350^2+AT1G09540+AT1G09540^2+AT2G18060+AT2G18060^2 +AT1G03040+ AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040+AT5G67420^2+AT1G25550^2+AT1G68670^2+AT1G66140^2 +AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';


% modelspec5 = 'AT1G08090_1~LBD37*HHO3*HHO2*ZFP4*BZIP53*NAP57*AT4G12750*BZO2H1*XERICO*AT2G43500*MYB42*MYB61*VND1*AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';

modelspec5 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';

% modelspec5 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040';
% modelspec = 'AT1G08090_1~ AT5G67420 * AT5G67420^2 * AT1G25550 * AT1G25550^2 *AT1G68670*AT1G68670^2 *AT1G66140*AT1G66140^2 *AT3G62420*AT3G62420^2 *AT3G57150*AT3G57150^2 *AT4G12750*AT4G12750^2*AT4G02640^2*AT4G02640*AT2G04240*AT2G04240^2*AT2G43500*AT2G43500^2*AT4G12350*AT4G12350^2*AT1G09540*AT1G09540^2*AT2G18060*AT2G18060^2 *AT1G03040*AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550^2+AT1G68670^2+AT1G66140^2+AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550+AT1G68670+AT1G66140+AT3G62420+AT3G57150+AT4G12750+AT4G02640+AT2G04240+AT2G43500+AT4G12350+AT1G09540+AT2G18060+AT1G03040';

% mdl = fitlm(Expmt,'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040');
% mdl = fitlm(Expmt,modelspec);
mdl5 = fitlm(Expmt5,modelspec5);
eee = find(mdl5.Coefficients.Estimate)
name555 = mdl5.CoefficientNames(eee)
model50 = mdl5.Coefficients.Estimate(eee);
figure()
h5 = bar(mdl5.Coefficients.Estimate(eee));
set(h5,'facecolor',[0.8 0.8 0.9])
legend('Coefficient')
set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
set(gca,'xticklabel',[name555(1) name555(2) name555(3) name555(4) name555(5) name555(6) name555(7) name555(8) name555(9) name555(10)])
ylabel('Input Genes')
xlabel('Normalized Coefficient')
title('Quadratic Model Coefficients')

Experiment_555 = table(name555',model50);
Experiment_555.Properties.VariableNames = {'Gene' 'Coefficient'};
writetable(Experiment_555,'Experiment 555.xlsx')

% 
% plotSlice(mdl5)

%%
x6 = readtable('out6.csv')

input_gene6 = x6(:,(1:14))
output_gene6 = x6(:,15)
a6 = table2array(input_gene6(:,1))
b6= table2array(input_gene6(:,2))
c6= table2array(input_gene6(:,3))
d6= table2array(input_gene6(:,4))
e6= table2array(input_gene6(:,5))
f6= table2array(input_gene6(:,6))
g6= table2array(input_gene6(:,7))
h6= table2array(input_gene6(:,8))
i6= table2array(input_gene6(:,9))
j6= table2array(input_gene6(:,10))
l6= table2array(input_gene6(:,11))
m6= table2array(input_gene6(:,12))
n6= table2array(input_gene6(:,13))
o6 = table2array(input_gene6(:,14))

input_gene6 = [a6 b6 c6 d6 e6 f6 g6 h6 i6 j6 l6 m6 n6 o6];
input_gene6 = normalize(input_gene6);
output_gene6 = table2array(output_gene6);
% input_gene = table2array(input_gene);
% output_gene = table2array(output_gene);
runorder = randperm(10);
Expmt6 = table(runorder', input_gene6(:,1), input_gene6(:,2), input_gene6(:,3), input_gene6(:,4), input_gene6(:,5), input_gene6(:,6),input_gene6(:,7), input_gene6(:,8),...
    input_gene6(:,9),input_gene6(:,10),input_gene6(:,11),input_gene6(:,12),input_gene6(:,13),input_gene6(:,14), output_gene6(:,1),'VariableNames',{'RunNumber','LBD37','HHO3','HHO2','ZFP4','BZIP53','NAP57','AT4G12750','BZO2H1','XERICO','AT2G43500','MYB42','MYB61','VND1','AT1G03040','AT1G08090_1'});
% Expmt = table(runorder', input_gene(:,1), input_gene(:,2), input_gene(:,3), input_gene(:,4), input_gene(:,5), input_gene(:,6),input_gene(:,7), input_gene(:,8),...
%     input_gene(:,9),input_gene(:,10),input_gene(:,11),input_gene(:,12),input_gene(:,13),input_gene(:,14), output_gene(:,1),'VariableNames',{'RunNumber','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'});
% mdl = fitlm(Expmt,'o~a*b*c*d*e*f*g*h*i*j*k*l*m*n');

% modelspec6 = 'AT1G08090_1~LBD37*HHO3*HHO2*ZFP4*BZIP53*NAP57*AT4G12750*BZO2H1*XERICO*AT2G43500*MYB42*MYB61*VND1*AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';
modelspec6 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040+LBD37^2+HHO3^2+HHO2^2+ZFP4^2+BZIP53^2+NAP57^2+AT4G12750^2+BZO2H1^2+XERICO^2+AT2G43500^2+MYB42^2+MYB61^2+VND1^2+AT1G03040^2';
% modelspec6 = 'AT1G08090_1~LBD37+HHO3+HHO2+ZFP4+BZIP53+NAP57+AT4G12750+BZO2H1+XERICO+AT2G43500+MYB42+MYB61+VND1+AT1G03040';

% modelspec = 'AT1G08090_1~ AT5G67420 * AT5G67420^2 * AT1G25550 * AT1G25550^2 *AT1G68670*AT1G68670^2 *AT1G66140*AT1G66140^2 *AT3G62420*AT3G62420^2 *AT3G57150*AT3G57150^2 *AT4G12750*AT4G12750^2*AT4G02640^2*AT4G02640*AT2G04240*AT2G04240^2*AT2G43500*AT2G43500^2*AT4G12350*AT4G12350^2*AT1G09540*AT1G09540^2*AT2G18060*AT2G18060^2 *AT1G03040*AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550^2+AT1G68670^2+AT1G66140^2+AT3G62420^2+AT3G57150^2+AT4G12750^2+AT4G02640^2+AT2G04240^2+AT2G43500^2+AT4G12350^2+AT1G09540^2+AT2G18060^2+AT1G03040^2';
% modelspec = 'AT1G08090_1~AT5G67420+AT1G25550+AT1G68670+AT1G66140+AT3G62420+AT3G57150+AT4G12750+AT4G02640+AT2G04240+AT2G43500+AT4G12350+AT1G09540+AT2G18060+AT1G03040';

% mdl = fitlm(Expmt,'AT1G08090_1~AT5G67420*AT1G25550*AT1G68670*AT1G66140*AT3G62420*AT3G57150*AT4G12750*AT4G02640*AT2G04240*AT2G43500*AT4G12350*AT1G09540*AT2G18060*AT1G03040');
% mdl = fitlm(Expmt,modelspec);
mdl6 = fitlm(Expmt6,modelspec6);
fff = find(mdl6.Coefficients.Estimate)
model60 = mdl6.Coefficients.Estimate(fff)';
name666 = mdl6.CoefficientNames(fff)
figure()
h6 = bar(mdl6.Coefficients.Estimate(fff));
set(h6,'facecolor',[0.8 0.8 0.9])
legend('Coefficient')
set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
set(gca,'xticklabel',[name666(1) name666(2) name666(3) name666(4) name666(5) name666(6) name666(7) name666(8) name666(9) name666(10)])
ylabel('Input Genes')
xlabel('Normalized Coefficient')
title('Quadratic Model Coefficients')

Experiment_66 = table(name666',model60');
Experiment_66.Properties.VariableNames = {'Gene' 'Coefficient'};
writetable(Experiment_66,'Experiment 66.xlsx')

% plotSlice(mdl6)

%%
%this section is to combine all coefficients of all 6 experiments, and get
%the mean and standard deviation of the coefficients, to generalize the
%sensitivity analysis 


% model1 = mdl.Coefficients.Estimate;
% model2 = mdl2.Coefficients.Estimate;
% model3 = mdl3.Coefficients.Estimate;
% model4 = mdl4.Coefficients.Estimate;
% model5 = mdl5.Coefficients.Estimate;
% model6 = mdl6.Coefficients.Estimate;
% 
% 
% aa = find(model1)
% bb = find(model2)
% cc = find(model3)
% dd = find(model4)
% ee = find(model5)
% ff = find(model6)
% 
% 
% name1 = mdl.CoefficientNames;
% name2 = mdl2.CoefficientNames;
% name3 = mdl3.CoefficientNames;
% name4 = mdl4.CoefficientNames;
% name5 = mdl5.CoefficientNames;
% name6 = mdl6.CoefficientNames;
% 
% name1 = name1(aa);
% name2 = name2(bb);
% name3 = name3(cc);
% name4 = name4(dd);
% name5 = name5(ee);
% name6 = name6(ff);
% 
% model11=model1(aa);
% model22 = model2(bb);
% model33= model3(cc);
% model44=model4(dd);
% model55=model5(ee);
% model66=model6(ff);
% coeff1 = [model1(aa(1)) model2(bb(1)) model3(cc(1)) model4(dd(1)) model5(ee(1)) model6(ff(1))]; 
% coeff2 = [model1(aa(2)) model2(bb(2)) model3(cc(2)) model4(dd(2)) model5(ee(2)) model6(ff(2))];
% coeff3 = [model1(aa(3)) model2(bb(3)) model3(cc(3)) model4(dd(3)) model5(ee(3)) model6(ff(3))];
% coeff4 = [model1(aa(4)) model2(bb(4)) model3(cc(4)) model4(dd(4)) model5(ee(4)) model6(ff(4))];
% coeff5 = [model1(aa(5)) model2(bb(5)) model3(cc(5)) model4(dd(5)) model5(ee(5)) model6(ff(5))];
% coeff6 = [model1(aa(6)) model2(bb(6)) model3(cc(6)) model4(dd(6)) model5(ee(6)) model6(ff(6))];
% coeff7 = [model1(aa(7)) model2(bb(7)) model3(cc(7)) model4(dd(7)) model5(ee(7)) model6(ff(7))];
% coeff8 = [model1(aa(8)) model2(bb(8)) model3(cc(8)) model4(dd(8)) model5(ee(8)) model6(ff(8))];
% coeff9 = [model1(aa(9)) model2(bb(9)) model3(cc(9)) model4(dd(9)) model5(ee(9)) model6(ff(9))];
% coeff10 = [model1(aa(10)) model2(bb(10)) model3(cc(10)) model4(dd(10)) model5(ee(10)) model6(ff(10))];
% 
% 
% 
% % coeff1 = [model1(aa(1)) model1(aa(2)) model1(aa(3)) model1(aa(4)) model1(aa(5)) model1(aa(6)) model1(aa(7)) model1(aa(8)) model1(aa(9)) model1(aa(10))]';
% % coeff2 = [model2(bb(1)) model2(bb(2)) model2(bb(3)) model2(bb(4)) model2(bb(5)) model2(bb(6)) model2(bb(7)) model2(bb(8)) model2(bb(9)) model2(bb(10))]';
% % coeff3 = [model3(cc(1)) model3(cc(2)) model3(cc(3)) model3(cc(4)) model3(cc(5)) model3(cc(6)) model3(cc(7)) model3(cc(8)) model3(cc(9)) model3(cc(10))]';
% % coeff4 = [model4(dd(1)) model4(dd(2)) model4(dd(3)) model4(dd(4)) model4(dd(5)) model4(dd(6)) model4(dd(7)) model4(dd(8)) model4(dd(9)) model4(dd(10))]';
% % coeff5 = [model5(ee(1)) model5(ee(2)) model5(ee(3)) model5(ee(4)) model5(ee(5)) model5(ee(6)) model5(ee(7)) model5(ee(8)) model5(ee(9)) model5(ee(10))]';
% % coeff6 = [model6(ff(1)) model6(ff(2)) model6(ff(3)) model6(ff(4)) model6(ff(5)) model6(ff(6)) model6(ff(7)) model6(ff(8)) model6(ff(9)) model6(ff(10))]';
% 
% 
% 
% avg1 = mean(coeff1);
% avg2 = mean(coeff2);
% avg3 = mean(coeff3);
% avg4 = mean(coeff4);
% avg5 = mean(coeff5);
% avg6 = mean(coeff6);
% avg7 = mean(coeff7);
% avg8 = mean(coeff8);
% avg9 = mean(coeff9);
% avg10 = mean(coeff10);
% 
% st1 = std(coeff1);
% st2 = std(coeff2);
% st3 = std(coeff3);
% st4 = std(coeff4);
% st5 = std(coeff5);
% st6 = std(coeff6);
% st7 = std(coeff7);
% st8 = std(coeff8);
% st9 = std(coeff9);
% st10 = std(coeff10);
% 
% err_bar = [st1 st2 st3 st4 st5 st6 st7 st8 st9 st10];
% global_average = [avg1 avg2 avg3 avg4 avg5 avg6 avg7 avg8 avg9 avg10];
% 
% % Expmt60 = table(global_average(1), global_average(2), global_average(3), global_average(4), global_average(5),global_average(6),global_average(7), global_average(8),...
% %     global_average(9),global_average(10),'VariableNames',{'MYB42','HHO3_HHO2','BZO2H1_AT2G43500','VND1_AT1G03040','AT4G12750_BZO2H1_AT1G03040','BZIP53_HHO3_NAP57_AT4G12750','BZIP53/AT4G12750/BZO2H1/XERICO/MYB61','XERICO_AT4G12750_MYB42_MYB61_VND1','LBD37_HHO3_HHO2_ZFP4_ZFP4_AT2G43500_MYB42_MYB61','LBD37_HHO3_HHO2_ZFP4_NAP57_BZIP53_MYB42_MYB61_AT1G03040'});
% % 
% % figure()
% % hold on
% % h6 = bar(1:10,global_average);
% % errorbar(1:10,global_average,err_bar,'.')
% % set(h6,'facecolor',[0.8 0.8 0.9])
% % legend('Coefficient')
% % set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
% % set(gca,'xticklabel',[name1((1)) name1((2)) name1((3)) name1((4)) name1((5)) name1((6)) name1((7)) name1((8)) name1((9)) name1((10))], 'fontsize', 8)
% % ylabel('Normalized Coefficient')
% % xlabel('Gene')
% % title('Quadratic Model Coefficients')
% % hold on
% 
% 
% % Expmt60 = table(global_average(1), global_average(2), global_average(3), global_average(4), global_average(5),global_average(6),global_average(7), global_average(8),...
% %     global_average(9), global_average(10), 'VariableNames',{name1((1)), name1((2)), name1((3)), name1((4)), name1((5)), name1((6)), name1((7)), name1((8)), name1((9)), name1((10))});
% 
% Expmt61 = table(name1',global_average');
% Expmt61.Properties.VariableNames = {'Gene' 'Coefficient'};
% writetable(Expmt61,'gene coefficients.xlsx')