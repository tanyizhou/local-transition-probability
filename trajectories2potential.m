%clc; clear; close all
%24/07/2017

%trajectories_DataName = [foldname '\Track_' DataName '.csv'];
min_t = 0;
max_t = 1;
% min_pixel = -40;
% max_pixel = 40;
%bin_size = 1;
%N = (max_pixel - min_pixel)/bin_size + 1;
timestep = 1;
%
%%
%import_trajectories
%alldata = xlsread(trajectories_DataName);
alldata = histogram;
index = -1;
trajectories = cell(size(unique(alldata(:,1))));
[C,ia,ic] =unique(alldata(:,1));
alldata(:,3) = round(alldata(:,3)/bin_size);

for i = 1:size(trajectories,1) - 1
    trajectories{i} = alldata(ia(i):ia(i+1)-1,:);
end
trajectories{end} = alldata(ia(end):end,:);

%%
%count_local_left_exits
min_bin = round(min_pixel/bin_size);
max_bin = round(max_pixel/bin_size);
N = max_bin - min_bin +1;
count = zeros(N,1);
left_exits = zeros(N,1);
total_exits = zeros(N,1);

for i = 1:size(trajectories,1)
    Trajectories = trajectories{i};
    start = Trajectories(1,3);
    start_in_range = (min_bin<start)&&(start<max_bin);
    for j = 1:size(Trajectories,1) - 1
        if start_in_range
            count(start - min_bin + 1) = count(start - min_bin + 1) + 1;
        end
        if Trajectories(j,3) ~= start
            if start_in_range
                count(start - min_bin + 1) = count(start - min_bin + 1) + 1;
                if Trajectories(j,3) < start
                    left_exits(start - min_bin + 1) = left_exits(start - min_bin + 1) + 1;
                end
                total_exits(start - min_bin + 1) = total_exits(start - min_bin + 1) + 1;
            end
            start = Trajectories(j,3);
            start_in_range = (min_bin<start)&&(start<max_bin);
        end
    end   
end

Dict = [left_exits, count, total_exits]; %total_counts is count
%DataName = [DataName '.mat'];
%save(DataName,'Dict','bin_size','diffusivity','timestep','N','min_pixel','max_pixel');

%%
splitting_prob_model = {
  'data{'
  'int N;'
  'int total_counts[N];'
  'int total_exits[N];'
  'int left_exits[N];'
  'real diffusivity;'
  'real timestep;'
  'real bin_size;'
  '}'
  'parameters {'
  'vector<lower=-1, upper=1>[N] F; // kB T = 1'
'}'

'transformed parameters {'
  'vector<lower=0, upper=1>[N] left_probability;'
  'vector<lower=-1.2/bin_size, upper=1.2/bin_size>[N] kappa;'
  'vector[N] U;'

  '// U[-1] = 0 by definition'
  'U[1] = -F[1] * bin_size;'
  'for (n in 2:N) {'
    'U[n] = U[n-1] - F[n] * bin_size;'
  '}'
  'for (n in 2:(N-1)) {'
    'kappa[n] = 0.5*(F[n+1] - F[n-1])/bin_size;'
  '}'
  'kappa[1] = kappa[2];'
  'kappa[N] = kappa[N-1];'
  'for (n in 1:N) {'
    'if(fabs(kappa[n]) > 0.05*fabs(F[n])) {'
      'left_probability[n] = (erfc((pow(2,-0.5)*(bin_size + 2*(1 - exp(-(diffusivity*timestep*kappa[n])))*F[n]*pow(kappa[n],-1))*pow((1 - exp(-2*diffusivity*timestep*kappa[n]))*pow(kappa[n],-1),-0.5))/2.)*'
        'pow(0.5 + erfc((pow(2,-0.5)*(bin_size + 2*(1 - exp(-(diffusivity*timestep*kappa[n])))*F[n]*pow(kappa[n],-1))*pow((1 - exp(-2*diffusivity*timestep*kappa[n]))*pow(kappa[n],-1),-0.5))/2.)/2. -'
        '(erf((diffusivity*pow(pow(diffusivity,-2)*pow(kappa[n],-1)*pow(2*F[n] + exp(diffusivity*timestep*kappa[n])*(-2*F[n] + bin_size*kappa[n]),2)*'
        '(-1 + cosh(diffusivity*timestep*kappa[n])*pow(sinh(diffusivity*timestep*kappa[n]),-1)),0.5))/4.)*(bin_size/2. + (-1 + exp(-(diffusivity*timestep*kappa[n])))*F[n]*pow(kappa[n],-1))*'
       ' pow(pow(bin_size/2. + (-1 + exp(-(diffusivity*timestep*kappa[n])))*F[n]*pow(kappa[n],-1),2),-0.5))/2.,-1))/2.;'
   ' } else {'
      'left_probability[n] = (exp(-(pow(bin_size,2)*pow(diffusivity,-1)*pow(timestep,-1))/8. - (diffusivity*timestep*pow(F[n],2))/2.)*pow(pi(),-0.5)*'
        '(4*erfc(((bin_size + 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.)*'
        '(erfc(((bin_size - 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.) + erfc(((bin_size + 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.))*'
        'exp((pow(bin_size,2)*pow(diffusivity,-1)*pow(timestep,-1))/8. + (diffusivity*timestep*pow(F[n],2))/2.)*pow(pi(),0.5) -'
        'bin_size*erfc(((bin_size - 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.)*exp((pow(diffusivity,-1)*pow(timestep,-1)*pow(bin_size - 2*diffusivity*timestep*F[n],2))/16.)*kappa[n]*'
        'pow(diffusivity*timestep,0.5) + bin_size*erfc(((bin_size + 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.)*'
        'exp((pow(diffusivity,-1)*pow(timestep,-1)*pow(bin_size + 2*diffusivity*timestep*F[n],2))/16.)*kappa[n]*pow(diffusivity*timestep,0.5))*'
        'pow(erfc(((bin_size - 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.) + erfc(((bin_size + 2*diffusivity*timestep*F[n])*pow(diffusivity*timestep,-0.5))/4.),-2))/4.;'
   ' }'
'  }'
'}'
'model {'
  'for (n in 1:N) {'
    'left_exits[n] ~ binomial(total_exits[n], left_probability[n]);'
  '}'
'}'
};


splitting_prob_dat = struct('N', N,...
                                   'total_counts', Dict(:,2),...
                                   'total_exits', Dict(:,3),...
                                   'left_exits', Dict(:,1),...
                                   'diffusivity', diffusivity,...
                                   'timestep',timestep,...
                                   'bin_size', bin_size);
                               
fit = stan('model_code',splitting_prob_model,'data',splitting_prob_dat, 'chains',1,'file', 'sp_model.stan');
fit.verbose = true;
fit.block();
%%
Temp = fit.print();

%dlmcell(DataName,Temp);


x = (min_pixel:bin_size:max_pixel);
U = fit.extract('permuted',true).U;
Umean = mean(U)';
Ustd = std(U)';
%p = mfilename('fullpath') ;
figure('Name',p);
shadedErrorBar(x, Umean, Ustd,{'k','LineWidth',1},1); hold on;
%maxUsp = Umean;
%Potential = Potentail + 6;
Potential = Potential + (Umean(round(size(Umean,1)/2))-Potential(round(size(Potential,1)/2)));
plot(Position, Potential, 'b', 'LineWidth', 2);
title(['Min=' num2str(min_pixel) ',Max=' num2str(max_pixel) ',size=' num2str(bin_size)]);
xlabel(DataName,'Interpreter', 'none');
dlmwrite([foldname '\SP_' DataName '.txt'],[x', Umean, Ustd]);

%print('-clipboard','-dmeta')