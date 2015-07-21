FEALevelSetWrapperMF_v3(5);

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = [0 0 3.25 3];
fig.PaperPositionMode = 'manual';


print('fig1_5_points','-dpng')
FEALevelSetWrapperMF_v3(10);
print('fig2_10_points','-dpng', '-r600')