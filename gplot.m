%gplot plots all the measures calculated in GrowMeasures
%function gplot(co1)

col1=['.' c1]

subplot (4,2,1)
plot(s.cost,mean(s.arand,1),'.k');
hold on
plot(s.cost,mean(s.a,1),col1);
ylim([-0.5 1])
xlabel('Cost');
ylabel('a');
title 'Assortativity'

subplot (4,2,2)
plot(s.cost,mean(s.Mrand,1),'.k');
hold on
plot(s.cost,mean(s.M,1),col1)
ylim([-0.5 1])
xlabel('Cost');
ylabel('M');
title 'Modularity'

subplot (4,2,3)
plot(s.cost,mean(s.Crand,1),'.k')
hold on
plot(s.cost,mean(s.C,1),col1)
xlabel('Cost');
ylabel('C');
title 'Clustering Coef'

subplot (4,2,4)
plot(s.cost,mean(s.Lrand,1),'.k')
hold on
plot(s.cost,mean(s.L,1),col1)
xlabel('Cost');
ylabel('L');
title 'Path Length'

subplot (4,2,5)
plot(s.cost,ones(length(s.cost),1),'.k')
hold on
plot(s.cost,mean(s.Sigma,1),col1)
ylim([0 10])
xlabel('Cost');
ylabel('Sigma');
title 'Small-World Coef'

subplot (4,2,6)
plot(s.cost,mean(s.Erand,1),'.k');
hold on
plot(s.cost,mean(s.E,1),col1)
xlabel('Cost');
ylabel('E');
title 'Efficiency'

subplot (4,2,7)
plot(s.cost,mean(s.CErand,1),'.k');
hold on
plot(s.cost,mean(s.CE,1),col1)
ylim([0 0.6])
xlabel('Cost');
ylabel('CE');
title 'Cost Efficiency'

