%% fitsensor: use fitting function to remove baseline shift
% % Inputs:
% data1: coG
% data2: coTd
% degree: fiting degree, use 1 for least square fitting

% % Outputs:
% p: fitted curve parameters 
% Linear model Poly1:
% f(x) = p(1,1)*x + p(1,2)
% p1: fitted curve parameters for coG
% p2: fitted curve parameters for coTd

function [p1,p2,data1_fit,data2_fit,ratio]= fitsensor (data1,data2,degree) 
    t = 1:1:length(data1);                  % generate a time scale
    p1 = polyfit(t,data1,degree);           % polyfit data1 to time
    data1_p1 = polyval(p1,t);               % use polyfit model, calculate the fitted values
    data1_fit = data1-data1_p1+mean(data1);             % subtract fitting lines from raw data,
                                            % add back mean of raw data for the y shift
      
    p2 = polyfit(t,data2,degree);
    data2_p2 = polyval(p2,t);
    data2_fit = data2-data2_p2+mean(data2);
    
    ratio = data1_fit./data2_fit;           % normalize data1(coG) to data2(coTd) for 
                                            % the actual signal changes
                                           
   figure('Position', [100 100 900 600])    % plot everything for double check
    subplot(3,1,1)
    plot(t,data1,'k')
    hold on
    plot(t,data1_p1,'b')
    hold on
    plot(t,data1_fit,'r')
    legend('raw','polyfit','fitted data','Location','eastoutside')
    legend('boxoff')
    title(['coG' '    y=' num2str(p1(1)) '*x+' num2str(p1(2))])

    subplot(3,1,2)
    plot(t,data2,'k')
    hold on
    plot(t,data2_p2,'b')
    hold on
    plot(t,data2_fit,'r')
    legend('raw','polyfit','fitted data','Location','eastoutside')
    legend('boxoff')
    title(['coTd' '    y=' num2str(p2(1)) '*x+' num2str(p2(2))])

    subplot(3,1,3)
    plot(t,data1./data2,'k')
    hold on
    plot(t,ratio,'r')
    legend('raw','corrected','Location','eastoutside')
    legend('boxoff')
    title('coG coTd ratio')
    
    
end
