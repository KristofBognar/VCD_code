function print_setting(size,save_fig,fig_name)
% size, size of the figure: size = 1, 1/2, 1/4, 0, or 'XL'
% save_fig, save the figure or not? : 1 = yes; 0 = no
% fig_name, name of the figure: any string you want, no extension!

if nargin == 2
    fig_name = 'test_output'; % if no figure name from input, then name it as "test_output"
elseif nargin == 1
    fig_name = 'test_output';% if no specify save figure or not, default is not to save!
    save_fig = 0;
end

switch size % define the size of the figure
    case 1
        W=19.0;H=23; %1 page figure
    case 1/2
        W=19.0;H=11.05; %1/2 figure
    case 1/4
        W=11.5;H=9.5; %1/4 figure
    case 0
        W=16.8 %cm
        H=16.8*3/4 %cm
    case 'XL'
        W = 25; H = 25;
end

set(gcf,'Units','centimeters');
set(gcf,'Position',[3 4 W H]);
set(gcf,'paperunits','centimeters');
set(gcf,'papersize',[W H]);
set(gcf,'paperposition',[3 4 W H]);
set(gca,'Fontname','Times New Roman','FontSize',14);
set(gcf,'PaperPositionMode', 'auto');
% print(gcf,'-depsc2','-loose','-r1200',['SOAZ.eps']);
% print(gcf,'-dpdf','-loose','-r1200',['2.jpg']);

if save_fig == 1 % define the output format
    saveas(gca,[fig_name, '.jpg']);
    saveas(gca,[fig_name, '.fig']);
    %saveas(gca,[fig_name, '.pdf']);
end