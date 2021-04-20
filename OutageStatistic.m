% compare GLINS .coor results with GINS reference result
% 10 or more columns of data should be input in format:
% sow lat lon height v_n v_e v_d roll pitch yaw ...
% unit: second, degree, m, m/s

%%
clear all;
close all;
clc;

%% input files: ground truth - ref_files, tested results - cmp_files, outages in GINS format - outage_files
ref_files = {'ref_file1.txt';
    'ref_file2.txt'};
cmp_files = {'cmp_file1.coor';
    'cmp_file2.coor'};
outage_files = {'outage1.outage';
    'outage2.outage'};

% variable to store results
MAX_POS_ERRS = [];
MAX_PLANE_ERRORS = [];
PLANE_DISPLACEMENTS = [];
MAX_ATT_ERRS = [];

% several results are supported
for i = 1:length(outage_files)
    skip_outage = [];
    
    cmp_data = importdata(cmp_files{i});
    ref_data = importdata(ref_files{i});
    
    % get the intersection part
    cmp_data = cmp_data(cmp_data(:,1) < ref_data(end,1),:);
    cmp_data = cmp_data(cmp_data(:,1) > ref_data(1,1),:);
    ref_id = find(ref_data(:,1) > cmp_data(1,1));
    ref_data = ref_data(ref_id,:);
    
    outages = importdata(outage_files{i});
    outages(:,3) = outages(:,2)+outages(:,3);
    
    % remove duplicate data and interpolate
    [~, index] = unique(ref_data(:, 1));
    ref_data = ref_data(index, :);
    ref_data(:,10) = smooth_angle(ref_data(:,10));
    
    [~, index] = unique(cmp_data(:, 1));
    cmp_data = cmp_data(index, :);
    
    tmp = zeros(size(cmp_data,1),10);
    tmp(:,1) = cmp_data(:,1);
    time = cmp_data(:,1);
    for j = 2:size(ref_data, 2)
        tmp(:, j) = interp1(ref_data(:, 1), ref_data(:, j), tmp(:,1));
    end
    ref_data = tmp;
    
    % coordinate transform
    a = 6378137.0;
    b = 6356752.3141;
    e2 = 0.00669437999013;
    
    ref_data(:, 2:3) = ref_data(:, 2:3) * pi / 180.0;
    tmp = 1-e2*sin(ref_data(:, 2)).^2;
    Rm1 = a * (1-e2) ./ tmp.^(3/2);
    Rn1 = a ./ sqrt(tmp);
    
    cmp_data(:, 2:3) = cmp_data(:, 2:3) * pi / 180.0;
    
    % statistic of errors
    pos_error = zeros(length(ref_data),3);
    pos_error(:,1) = (Rm1 + ref_data(:,4)).*(ref_data(:, 2) - cmp_data(:, 2)); % N error
    pos_error(:,2) = (Rn1 + ref_data(:,4)).*(ref_data(:, 3) - cmp_data(:, 3)).*cos(ref_data(:, 2)); % E error
    pos_error(:,3) = ref_data(:, 4) - cmp_data(:, 4); % D error
    att_error = ref_data(:, 8:10) - cmp_data(:, 8:10);
    att_error(:,end) = wrap_angle(att_error(:,end));
    
    % get mileage
    dif_phi = diff(ref_data(:,2));
    dif_lambda = diff(ref_data(:,3));
    dpos = zeros(length(dif_phi), 3);
    dpos(:,1) = (Rm1(1:end-1) + ref_data(1:end-1,4)).* dif_phi;
    dpos(:,2) = (Rn1(1:end-1) + ref_data(1:end-1,4)).* dif_lambda.* cos(ref_data(1:end-1,2));
    dpos(:,3) = diff(ref_data(:,4));
    
    for j = 1:size(outages,1)
        idx = find(time >= outages(j,2) & time <= outages(j,3));
        
        if isempty(idx)
            skip_outage = [skip_outage, j];
            continue;
        end
        
        if idx(1) > 1
            idx = idx-1;
        else
            idx = idx(2:end);
        end            
        
        pos_err_set = pos_error(idx,:);
        MAX_POS_ERRS = [MAX_POS_ERRS; max(abs(pos_err_set),[],1)];
        
        plane_err_set = sqrt(pos_err_set(:,1).^2+pos_err_set(:,2).^2);
        MAX_PLANE_ERRORS = [MAX_PLANE_ERRORS; max(plane_err_set)];
        
        plane_displacement = sum(sqrt(sum(dpos(idx,1:2).^2,2)));
        PLANE_DISPLACEMENTS = [PLANE_DISPLACEMENTS; plane_displacement];
        
        att_err_set = att_error(idx,:);
        MAX_ATT_ERRS = [MAX_ATT_ERRS; max(abs(att_err_set),[],1)];
    end
    
    max_y_pos = ceil(max(max(abs(MAX_POS_ERRS))));
    max_y_att = ceil(max(max(abs(MAX_ATT_ERRS))));
    
    % set xtick format
    nf = java.text.DecimalFormat;
    x_formatstring = '%6d';
    t_start = floor(cmp_data(1,1)/100)*100;
    
    % plot
    figure;
    set(gcf,'position',[150 10 660 550])
    subplot(211);
    hold on;
    plot(time-t_start, pos_error(:, 1), 'b.', 'MarkerSize', 10);
    plot(time-t_start, pos_error(:, 2), 'g.', 'MarkerSize', 10);
    plot(time-t_start, pos_error(:, 3), 'r.', 'MarkerSize', 10);
    
    for j = 1:size(outages,1)
        if ~isempty(find(skip_outage==j, 1))
            continue;
        end
        
        if mod(j,2) == 1
            xtick(uint32((j+1)/2)) = outages(j,2)-t_start;
%             xticklabel{uint32((j+1)/2)} = char(nf.format(xtick(uint32((j+1)/2))));
            xticklabel{uint32((j+1)/2)} = xtick(uint32((j+1)/2)); % for MDPI Sensors template
        end
        
        outage_start = outages(j, 2);
        outage_end = outages(j, 3);
        if outage_end > time(end)
            outage_end = time(end);
        end
        outage_range = [outage_start, outage_start, outage_end, outage_end];
        dy = [-max_y_pos,max_y_pos];
%         dy = [-5, 5];
        dy_range = [dy(1)+0.01, dy(2)-0.01, dy(2)-0.01, dy(1)+0.01];
        fill(outage_range-t_start, dy_range, [0.2, 0.2, 0.2], 'Edgecolor', 'none');
        alpha(0.2);
    end
    hold off;
    ylim(dy);
    l1 = legend('N', 'E', 'D', 'Location', 'eastoutside');
    set(l1,'Fontname', 'Palatino Linotype', 'box', 'off', 'FontSize', 15, 'FontWeight','bold');
    xlabel(['(sow - ', num2str(t_start), ') (s)'], 'FontName', 'Palatino Linotype', 'FontSize', 15, 'FontWeight','bold');
    ylabel('Position Error (m)', 'FontName', 'Palatino Linotype', 'FontSize', 15, 'FontWeight','bold');
    set(gca, 'xtick', xtick);
    set(gca, 'xticklabel', xticklabel, 'FontName', 'Palatino Linotype', 'FontSize', 15, 'FontWeight','bold');
    xticklabels(gca, strrep(xticklabels(gca),'-','–')); % from 'hyphen' to 'minus'
    yticklabels(gca, strrep(yticklabels(gca),'-','–')); % from 'hyphen' to 'minus'
    set(gca,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
    grid on;
    box on;
    
    subplot(212);
    hold on;
    plot(time-t_start, att_error(:, 1), 'b.', 'MarkerSize', 10);
    plot(time-t_start, att_error(:, 2), 'g.', 'MarkerSize', 10)
    plot(time-t_start, att_error(:, 3), 'r.', 'MarkerSize', 10);
    
    for j = 1:size(outages,1)
        if ~isempty(find(skip_outage==j, 1))
            continue;
        end
        
        outage_start = outages(j, 2);
        outage_end = outages(j, 3);
        if outage_end > time(end)
            outage_end = time(end);
        end
        outage_range = [outage_start, outage_start, outage_end, outage_end];
        dy = [-max_y_att,max_y_att];
%         dy = [-1, 1];
        dy_range = [dy(1)+0.01, dy(2)-0.01, dy(2)-0.01, dy(1)+0.01];
        fill(outage_range-t_start, dy_range, [0.2, 0.2, 0.2], 'Edgecolor', 'none');
        alpha(0.2);
    end
    hold off;
    ylim(dy);
    l2 = legend('R', 'P', 'Y', 'Location', 'eastoutside');
    set(l2,'Fontname', 'Palatino Linotype', 'box', 'off', 'FontSize', 15, 'FontWeight','bold');
    xlabel(['(sow - ', num2str(t_start), ') (s)'], 'FontName', 'Palatino Linotype', 'FontSize', 15, 'FontWeight','bold');
    ylabel('Attitude Error (\circ)', 'FontName', 'Palatino Linotype', 'FontSize', 15, 'FontWeight','bold', 'FontWeight','bold');
    set(gca, 'xtick', xtick);
    set(gca, 'xticklabel', xticklabel, 'FontName', 'Palatino Linotype', 'FontSize', 15, 'FontWeight','bold');
    xticklabels(gca, strrep(xticklabels(gca),'-','–')); % from 'hyphen' to 'minus'
    yticklabels(gca, strrep(yticklabels(gca),'-','–')); % from 'hyphen' to 'minus'
    set(gca,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
    grid on;
    box on;
end

figure;
b1 = bar(PLANE_DISPLACEMENTS, 'r', 'EdgeColor', 'None');
axis tight;
xlabel('Id of outages');
ylabel('Mileage during outages/m'); 
set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 15);

disp('pole statistic:');
str = sprintf('rms pos err in N, E, D:\n %.2f, %.2f, %.2f\n', rms(MAX_POS_ERRS,1));
disp(str);
str = sprintf('max pos err in N, E, D:\n %.2f, %.2f, %.2f\n', max(MAX_POS_ERRS,[],1));
disp(str);

str = sprintf('rms att err in R, P, Y:\n %.2f, %.2f, %.2f\n', rms(MAX_ATT_ERRS,1));
disp(str);
str = sprintf('max att err in R, P, Y:\n %.2f, %.2f, %.2f\n', max(MAX_ATT_ERRS,[],1));
disp(str);

str = sprintf('max pos err in plane: %.2f', max(MAX_PLANE_ERRORS));
disp(str);

% ref displacements is accuracy
str = sprintf('mean relative pos err in plane: %.3f%%', sum(MAX_PLANE_ERRORS)/sum(PLANE_DISPLACEMENTS)*100);
disp(str);