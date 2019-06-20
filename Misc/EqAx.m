%% Equalizes Axes

ehv=[xlim;ylim]';
xlim([min(ehv(:)) max(ehv(:))]);
ylim([min(ehv(:)) max(ehv(:))]);