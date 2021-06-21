function  c=Empirical_ccdf(x)
x = reshape(x,numel(x),1);%reshape x to [numel(x),1]
% select method (discrete or continuous)  
if isempty(setdiff(x,floor(x)))
    f_dattype = 'INTS';
elseif isreal(x)    
    f_dattype = 'REAL';
else
    f_dattype = 'UNKN';
end;


% estimate xmin and alpha, accordingly
switch f_dattype    
    case 'REAL'	
        n = length(x);
        c= [sort(x) (n:-1:1)'./n];
	case 'INTS'
        n = length(x);        
        q = unique(x);
        c = hist(x,q)'./n;        
        c = [[q; q(end)+1] 1-[0; cumsum(c)]];       
        c(c(:,2)<10^-10,:) = [];
    otherwise
        fprintf('Error: x must contain only reals or only integers.\n');
        return;
end

       