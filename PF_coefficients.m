function coef = PF_coefficients(input,bins)
%Applies a planar fit to force vbar = wbar = 0. u, v, and w must be contained
%in columns 1, 2, and 3 respectively.

numBins = length(bins);
if ~numBins; numBins = 1; end

for i = 1:numBins
    
    % sort data by direction and store in matrix 
    if i < length(bins)
        matrix = input(input(:,4)>bins(i) & input(:,4)<=bins(i+1),:);
    elseif isempty(bins)
        matrix = input;
    else
        matrix = input(input(:,4)>bins(i) | input(:,4)<=bins(1),:);
    end
    
    %Partition Data
    ubar = matrix(:,1);
    vbar = matrix(:,2);
    wbar = matrix(:,3);
    
    %Remove NaNs
    nanFlag = logical(isnan(ubar)+isnan(vbar)+isnan(wbar));
    ubar(nanFlag) = [];
    vbar(nanFlag) = [];
    wbar(nanFlag) = [];

    %Number of Records
    flen = length(ubar);
    
    %Find Values to Populate Matrix Equation
    su = sum(ubar);
    sv = sum(vbar);
    sw = sum(wbar);
    suv = sum(ubar.*vbar);
    suw = sum(ubar.*wbar);
    svw = sum(vbar.*wbar);
    su2 = sum(ubar.*ubar);
    sv2 = sum(vbar.*vbar);
    
    % [H](b) = (g)
    H = [flen su sv; su su2 suv; sv suv sv2];
    g = [sw suw svw]';
    
    %Solve for b matrix
    if i < length(bins)
        direction = sprintf('degrees_%g_to_%g',bins(i),bins(i+1));
        coef.(direction) = linsolve(H,g);
    elseif isempty(bins)
        direction = 'degrees_0_to_0';
        coef.(direction) = linsolve(H,g);
    else
        direction = sprintf('degrees_%g_to_%g',bins(i),bins(1));
        coef.(direction) = linsolve(H,g);
    end
end
end