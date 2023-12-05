%
% new_ages_ka = find_outliers_transect(ages_ka)
% new_ages_ka = find_outliers_transect(ages_ka,mask)
% new_ages_ka = find_outliers_transect(ages_ka,mask,strat_level)
% new_ages_ka = find_outliers_transect(ages_ka,mask,strat_level,exclude_ends)
%
% Finds outliers within a transect of exposure ages based on the
% stratigraphic relationship of the samples. The age of each sample (mean 
% % and 1 st.dev.) is compared to that of samples stratigraphically 
% above/below it. Each sample is assessed in turn, running from 
% top-to-bottom and then bottom-to-top within the transect. A sample is 
% considered a 'likely outlier' if it was relatively too old/young in 
% either direction, and a 'definite outlier' if relatively too old/young in 
% both directions.
%
% Note, currently set up for Be-10 exposure ages only.
%
% ages_ka is a required struct, containing ages and uncertainties 
% calculated using age_calc.m, elevations and positions, and a logical of 
% the nuclides that were measured. See iceTEA.
% 
% mask specifies which samples to include. This should be based on the 
% original input sample data. Leave empty ([]) to include all samples.
%
% strat_level specifies how many samples to consinder above and below each
% sample when assessing the stratigraphic relationship. This should be 2 or
% 3. The default is 3.
% 
% exclude_ends is used to specify whether to exclude the ages at either end
% of the transect (e.g. top and bottom samples) for the assessment - 
% 0=include, 1=exclude.
%
% Outputs a struct with the exposure age dataset (outliers removed), and 
% the likely and distinct outliers identified, if any.
%
%
% Written by Richard Selwyn Jones, Monash University
% richard.s.jones@monash.edu
% Compatible with the iceTEA tools suite.
%
%
%%

function new_ages_ka = find_outliers_transect(ages_ka,mask,strat_level,exclude_ends)

  % Check inputs
  if nargin > 4
      error('find_outliers_transect has wrong number of inputs!');
  end
  if (nargin < 2 || isempty(mask))
      mask = 1:length(ages_ka.sample_names);
  end
  if (nargin < 3 || isempty(strat_level))
      strat_level = 3;
  end
  if ~ismember(strat_level,[2,3])
      error('strat_level should be 2 or 3');
  end
  if (nargin < 4 || isempty(exclude_ends))
      exclude_ends = 0;
  end
  if ~ismember(exclude_ends,[0,1])
      error('exclude_ends should be binary [0 or 1]');
  end
  
  
  % Find and remove outliers for each nuclide
      
  if ages_ka.NN(1)
      
      sample_names = ages_ka.sample_names(mask);
      ages = ages_ka.Be10(mask,:);
      position = ages_ka.position(mask);
      if isnan(position)
          position = ages_ka.elevation(mask);
      end
      
      % Rank by relative position (largest to smallest)
      [srt_position,srt_pos_idx] = sort(position,'descend'); % Sort in descending order and get indices relative to original order
      [~, rnk_pos] = ismember(position,srt_position); % Rank in descending order, adjusting for any duplicate positions
      [uni_pos,~,rnk_pos_uni] = unique(srt_position,'stable'); % Rank in descending order for each unique position
      n=accumarray(rnk_pos_uni(:),1);  dup_vals=uni_pos(n~=1);  duplicate_pos_log = ismember(srt_position,dup_vals); % Get logical of duplicate positions
      
      srt_ages = ages(srt_pos_idx,:); % Sort ages based on position
      srt_names = sample_names(srt_pos_idx); % Sort sample names based on position
      
      srt_ages_rev = flipud(srt_ages);
      srt_position_rev = flipud(srt_position);
      srt_names_rev = flipud(srt_names);
      duplicate_pos_rev_log = fliplr(duplicate_pos_log);
      
      outlier_forward_log = false(1,length(srt_ages(:,1)));  outlier_reverse_log = outlier_forward_log;
            
      
      % Evaluate each sample position in forward stratigraphic order ('top-down')
      for a = 1:length(srt_ages(:,1))
          
          this_age = [srt_ages(a,1),srt_ages(a,3)]; % Get mean age with total uncertainty

          % Create new set of ages with any identified outliers removed
          clear temp_srt_ages
          temp_srt_ages = srt_ages(~outlier_forward_log,:);
          age_idx = find(this_age(1)==temp_srt_ages(:,1)); % Find age index within new set

          % Check age relative to next position (Level 1)
          if age_idx+1 <= numel(temp_srt_ages(:,1)) % Check there is a next position
              rel_out = age_rel(this_age,[temp_srt_ages(age_idx+1,1),temp_srt_ages(age_idx+1,3)]);

              if rel_out == 1 && age_idx+2 <= numel(temp_srt_ages(:,1)) && strat_level ~= 1
                  % If older than the next position, check age
                  % relative to the position after that (Level 2)
                  rel_out2 = age_rel(this_age,[temp_srt_ages(age_idx+2,1),temp_srt_ages(age_idx+2,3)]);

                  if rel_out2 == 3 && age_idx+3 <= numel(temp_srt_ages(:,1)) && strat_level == 3
                      % If now younger, check age relative to the position 
                      % after that (Level 3)
                      rel_out3 = age_rel(this_age,[temp_srt_ages(age_idx+3,1),temp_srt_ages(age_idx+3,3)]);

                      if rel_out3 == 3
                          this_outlier = 1; % If still younger, call it an outlier
                      else
                          this_outlier = 0; % If now older or consistent, cannot call it an outlier
                      end
                  else
                      this_outlier = 0; % If still older or consistent, cannot call it an outlier
                  end

              elseif rel_out == 3 && age_idx+2 <= numel(temp_srt_ages(:,1)) && strat_level ~= 1
                  % If younger than the next position, check age relative 
                  % to the position after that (Level 2)
                  rel_out2 = age_rel(this_age,[temp_srt_ages(age_idx+2,1),temp_srt_ages(age_idx+2,3)]);

                  if rel_out2 == 3
                      this_outlier = 1; % If also younger, call it an outlier
                  else
                      this_outlier = 0; % If older or consistent, do not call it an outlier
                  end
              else
                  this_outlier = 0; % If consistent with the next position, cannot call it an outlier
              end

          else
              if exclude_ends == 0
                  this_outlier = 1; % As the last position, unable to determine whether it is good (call it an outlier)
              elseif exclude_ends == 1
                  this_outlier = 0; % As the last position, unable to determine whether it is an outlier (call it good here)
              end
          end

          outlier_forward_log(a) = this_outlier;

      end
      
      
      % Evaluate each sample position in reverse stratigraphic order ('bottom-up')
      for b = 1:length(srt_ages_rev(:,1))

          this_age = [srt_ages_rev(b,1),srt_ages_rev(b,3)]; % Get mean age with total uncertainty

          % Create new set of ages with any identified outliers removed
          clear temp_srt_ages_rev
          temp_srt_ages_rev = srt_ages_rev(~outlier_reverse_log,:);
          age_idx = find(this_age(1)==temp_srt_ages_rev(:,1)); % Find age index within new set

          % Check age relative to previous position (Level 1)
          if age_idx+1 <= numel(temp_srt_ages_rev(:,1)) % Check there is a previous position
              rel_out = age_rel(this_age,[temp_srt_ages_rev(age_idx+1,1),temp_srt_ages_rev(age_idx+1,3)]);

              if rel_out == 3 && age_idx+2 <= numel(temp_srt_ages_rev(:,1)) && strat_level ~= 1
                  % If younger than the previous position, check age
                  % relative to the position before that (Level 2)
                  rel_out2 = age_rel(this_age,[temp_srt_ages_rev(age_idx+2,1),temp_srt_ages_rev(age_idx+2,3)]);

                  if rel_out2 == 1 && age_idx+3 <= numel(temp_srt_ages_rev(:,1)) && strat_level == 3
                      % If now older, check age relative to the position 
                      % before that (Level 3)
                      rel_out3 = age_rel(this_age,[temp_srt_ages_rev(age_idx+3,1),temp_srt_ages_rev(age_idx+3,3)]);

                      if rel_out3 == 1
                          this_outlier = 1; % If still older, call it an outlier
                      else
                          this_outlier = 0; % If now younger or consistent, cannot call it an outlier
                      end
                  else
                      this_outlier = 0; % If still younger or consistent, cannot call it an outlier
                  end

              elseif rel_out == 1 && age_idx+2 <= numel(temp_srt_ages_rev(:,1)) && strat_level ~= 1
                  % If older than the previous position, check age relative
                  % to the position before that (Level 2)
                  rel_out2 = age_rel(this_age,[temp_srt_ages_rev(age_idx+2,1),temp_srt_ages_rev(age_idx+2,3)]);

                  if rel_out2 == 1
                      this_outlier = 1; % If also older, call it an outlier
                  else
                      this_outlier = 0; % If younger or consistent, do not call it an outlier
                  end
              else
                  this_outlier = 0; % If consistent with the previous position, cannot call it an outlier
              end

          else
              if exclude_ends == 0
                  this_outlier = 1; % As the last position, unable to determine whether it is good (call it an outlier)
              elseif exclude_ends == 1
                  this_outlier = 0; % As the last position, unable to determine whether it is an outlier (call it good here)
              end
          end

          outlier_reverse_log(b) = this_outlier;

      end
           
      
      % Identify outliers and distinct outliers
      outlier_distinct_log = all([outlier_forward_log;fliplr(outlier_reverse_log)]); % Outliers when run in forward and reverse order
      outlier_likely_log = any([outlier_forward_log;fliplr(outlier_reverse_log)]) & ~outlier_distinct_log; % Outliers when run in just forward or reverse order
      
      % Get info on identified outliers
      out_likely_idx = find(outlier_likely_log);
      outliers_likely = srt_ages(out_likely_idx,1);
      outlier_likely_names = srt_names(out_likely_idx);
      out_distinct_idx = find(outlier_distinct_log);
      outliers_distinct = srt_ages(out_distinct_idx,1);
      outlier_distinct_names = srt_names(out_distinct_idx);
            
      
      % Export ages and outliers
      srt_ages(out_likely_idx,:) = NaN; % Remove any outliers from existing exposure age data
      new_ages_ka.Be10 = srt_ages; % Export ages without outliers
      
      new_ages_ka.Be10_outliers = outliers_likely; % Export likely outliers
      new_ages_ka.Be10_outliers_distinct = outliers_distinct; % Export distinct outliers
      
      
      % Find indices and values of outliers
      if ~isempty(outliers_distinct)
          disp('Distinct outliers (two-direction identification):')
          for c = 1:length(outliers_distinct)
              name = cellstr(outlier_distinct_names(c));
              disp(['sample ' name{1} ' (mean of ' sprintf('%0.2f',outliers_distinct(c)) ' ka)'])
          end
      end
      if ~isempty(outliers_likely)
          disp('Likely outliers (one-direction identification):')
          for c = 1:length(outliers_likely)
              name = cellstr(outlier_likely_names(c));
              disp(['sample ' name{1} ' (mean of ' sprintf('%0.2f',outliers_likely(c)) ' ka)'])
          end
      end
  end
  
  
end


%%%%%%%%%%%%% Determine age relationship between two samples %%%%%%%%%%%%%%

function out = age_rel(sample1_age,sample2_age)
  
  age1_low = sample1_age(1)-sample1_age(2);
  age1_upp = sample1_age(1)+sample1_age(2);
  
  age2_low = sample2_age(1)-sample2_age(2);
  age2_upp = sample2_age(1)+sample2_age(2);
  
  if age1_low > age2_upp
      out = 1; % Age 1 is older than age 2
  elseif age1_upp < age2_low
      out = 3; % Age 1 is younger than age 2
  else %(age1_upp > age2_low) && (age1_low < age2_upp)
      out = 2; % Age 1 and age 2 are consistent (overlap)
  end
  
end
