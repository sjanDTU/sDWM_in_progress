%% Efficiency versus iterations plot

idxEnd = find(fitMaxTrack==0,1)-1;

figure
hold on
plot(fitMaxTrack(1:idxEnd))
plot(fitMinTrack(1:idxEnd))
plot(fitMeanTrack(1:idxEnd))
xlabel('Iterations (-)')
ylabel('Total farm power (W)')

legend('Max. Weed','Min Weed','Mean Weed')

%% Operation strategy

figure
bar(fliplr(1-derOpt)*100)
xlabel('Turbine No. (-)')
ylabel('Derating (%)')


