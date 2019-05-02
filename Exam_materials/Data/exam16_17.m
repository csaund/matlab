%% MATLAB exam answers 2016-17


%% Q1. How do the 4 response types (hits, misses, false alarms, correct rejections) 
%% distribute within each culture for each facial expression type? 
% For each of the 15 individual subjects in each culture, compute the number of 
% hits, misses, false alarms, and correct rejections for each of the four facial expression types (thinking, interested, bored, confused)
% separately. Display the results of each culture as a 15 (subjects) X 4 (response types) color-coded matrix 
% where red indicates a high number and blue indicates a low number. To compare data across cultures, 
% use a 1 X 2 (cultures) subplot. Include all necessary features (e.g., title, axes titles, colorbar). (3 POINTS)


clear all 

% set up variables from folder 
% go to folder 
cd('/Users/carolynsaund/glasgow/matlab/Exam_materials/Data')
od = pwd;

% extract all variables from folders 
end_folder_name = ' subjects';
folders = dir(['*' end_folder_name]);
nr_cultures = length(folders);

cultures = cell(1, nr_cultures);
for cults = 1:nr_cultures
    cultures{cults} = folders(cults).name(1:end - length(end_folder_name));
end


% *assuming that all subjects have the same number of trials*
cd([od '/' cultures{1} end_folder_name])
files = dir('data_sj_*.mat');
nr_subjects = length(files);

load(files(1).name) 
tmp = lower(unique(word_all));
social_messages = tmp([4 3 1 2]); % re-order to match coding of responses
nr_social_messages = length(social_messages);
nr_trials = length(resps);


% now prepare the data in a useable way 
% prepare matrix for all data 
all_data = zeros(nr_trials, 4, nr_subjects, nr_cultures);

for culture = 1:nr_cultures
    
    % change directory for each culture
    cd([od '/' cultures{culture} end_folder_name])
    
    % directory of files
    files = dir('*.mat');
    
    for subject = 1:nr_subjects
        
       fname = files(subject).name;
       load(fname)
       
       % resps = 1 - yes, 2 - no
       % sti_cat = 1 - thinking, 2 - interested, 3 - bored, and 4 - confused
       % word_all = word presented on each trial
       
       % first, convert all words to numbers
       word_numbers = zeros(nr_trials, 1);
       for trials = 1:nr_trials
           word_numbers(trials, 1) = find(strcmpi(word_all{trials}, social_messages));
       end
       
       % responses,  facial expression types,  words presented on each trial
       all_data(:, :, subject, culture) = [sti_cat, word_numbers, resps, model_num];   
    end 
end 


%%

% hits, misses, false alarms, correct rejections 
response_types = {'hits', 'miss', 'FA', 'CR'};
nr_response_types = length(response_types); 
             

% analyse for each facial expression type and culture
% prepare matrix for response types
response_patterns = zeros(nr_subjects, nr_response_types, nr_social_messages, nr_cultures);

for facial_expression_type = 1:nr_social_messages  
    for culture = 1:nr_cultures
        for subject = 1:nr_subjects
            
            % extract relevant trials
            loc_trials = all_data(:, 1, subject, culture) == facial_expression_type; % facial expression type
            word_numbers = all_data(loc_trials, 2, subject, culture); % words
            resps = all_data(loc_trials, 3, subject, culture);  % yes/no
            
            % find and compute number of each response type
            hits = sum(resps == 1 & word_numbers == facial_expression_type); % hits
            misses = sum(resps == 2 & word_numbers == facial_expression_type); % misses
            false_alarms = sum(resps == 1 & word_numbers  ~= facial_expression_type); % false alarms
            correct_reject = sum(resps == 2 & word_numbers ~= facial_expression_type); % correct rejection
            
            % store data
            response_patterns(subject, :, facial_expression_type, culture) = [hits, misses, false_alarms, correct_reject];
        end    
    end 
end


%% display figure

% find max and min to display data using same scale
max_val = max(response_patterns(:));
min_val = min(response_patterns(:));

colormap jet
for facial_expression_type = 1:nr_social_messages
    figure,
    
    for culture = 1:nr_cultures
        
        subplot(1, 2, culture)
        imagesc(response_patterns(:, :, facial_expression_type, culture), [min_val max_val])
        set(gca, 'XTick', 1:nr_response_types, 'XTickLabel', response_types)
        title(cultures{culture})
        
        if culture == 1
            ylabel('Participants')
            xlabel('Response types')
        else
            colorbar
        end
    end
    % suptitle(social_messages{facial_expression_type});
end

   
%% Q2. Compute d-prime for each subject and facial expression type. 
% Using the necessary proportions of hits and false alarms, compute d-prime
% for each individual subject in each culture for each facial expression type separately. 
% Plot your results as individual points using a 2 x 2 (4 facial expression types) with the 
% data of each culture shown side by side. 
% Include all necessary features (e.g., title, axes titles, colorbar). (3 POINTS)


% d-prime is computed based on the proportion of hits and false alarms
% as computed across match and mismatch trials
% therefore, hits and false alarms are *independent*

% CARO adding this here because I guess somehow they put dprime in the wes
% sub folder?
cd('/Users/carolynsaund/glasgow/matlab/Exam_materials/Data')


% convert any 0s or 1s to avoif Inf
response_patterns(response_patterns == 1) = .99999;
response_patterns(response_patterns == 0) = .00001;


% prepare matrix for data storage
dp = zeros(nr_subjects, nr_social_messages, nr_cultures);
hits_fa_proportion = zeros(nr_subjects, 2, nr_social_messages, nr_cultures);

for culture = 1:nr_cultures
    for facial_expression_type = 1:nr_social_messages
        for subject = 1:nr_subjects
            
            % extract hits and false alarms
            hits = response_patterns(subject, 1, facial_expression_type, culture);
            false_alarms = response_patterns(subject, 3, facial_expression_type, culture);
            
            % compute proportion based on number of matches (hits) and
            % mismatches (false alarms)
            
            % get trial information
            % sti_cat, word_numbers, resps
            loc_trial_info = all_data(:, 1, subject, culture) == facial_expression_type;
            
            % count number of matches
            nr_matches = sum(all_data(loc_trial_info, 2, subject, culture) == facial_expression_type);
            nr_mismatches = sum(loc_trial_info) - nr_matches;
            
            % compute proportions
            hits_proportion = hits ./ nr_matches;
            FA_proportion = false_alarms ./ nr_mismatches;
            
            % store for later 
            hits_fa_proportion(subject, :, facial_expression_type, culture) = [hits_proportion, FA_proportion];
            
            % compute d-prime with downloaded function
            dp(subject, facial_expression_type, culture) = dprime_simple(hits_proportion, FA_proportion);
        end
    end
end


%% display figure 

% get max and min values to display on the same scale
max_dp = max(dp(:));
min_dp = min(dp(:));

figure, 
for facial_expression_type = 1:nr_social_messages
    for culture = 1:nr_cultures
        
        subplot(2, 2, facial_expression_type)
        
        data = dp(:, facial_expression_type, culture);
        y = ones(size(data)) .* culture;
        plot(y, data, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
        hold on 
        set(gca, 'XTick', 1:nr_cultures, 'XTickLabel', cultures)
        axis([0 nr_cultures + 1 min_dp max_dp + 1])
        ylabel('d-prime')
    end
    title(social_messages(facial_expression_type))
end

whitebg([0 0 0])


%% Q3. Do any subjects show a significantly higher or lower d-prime than the group (i.e., are potential outliers)? 
% For each culture and facial expression type separately, compute the Z-score of each subject's d-prime 
% and show whether it is significantly above, significantly below, or not significantly different 
% from the mean. Replot your data as in Q2 above using asterisks to indicate any 
% significant differences. Save this figure as a .tiff file in a new folder called Results Figures. (3 POINTS)


thresh = 1.96; % for two tailed 0.05
figure, 
for facial_expression_type = 1:nr_social_messages
    for culture = 1:nr_cultures
        
        data = dp(:, facial_expression_type, culture);
        % plot data as before
        
        subplot(2, 2, facial_expression_type)
        y = ones(size(data)) .* culture;
        plot(y, data, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
        hold on
        set(gca, 'XTick', 1:nr_cultures, 'XTickLabel', cultures)
        axis([0 nr_cultures + 1 min_dp max_dp + 1])
        ylabel('d-prime')
        hold on
        
        % Z score and significance
        Z_scores = (data - mean(data)) ./ std(data);
        signif_diffs = data(abs(Z_scores) >= thresh);
        nr_sig_diffs = length(signif_diffs);
        
        % plot asterisks
        for sig_diff = 1:nr_sig_diffs
            plot(culture, signif_diffs(sig_diff), '*-r')
            hold on 
        end 
    end
    title(social_messages(facial_expression_type))
end

% make new folder and save figure within it 
new_folder_name = 'Results Figures';
mkdir(new_folder_name)
cd([od '/' new_folder_name])
print('Distribution of d-prime values', '-dtiff')


%% Q4. Understanding confusions. 
% False alarms occur when subjects report the presence of a signal (here, a social message such as 'confused') when it is not present 
% (that is, when the facial expression presented does not match the label). To examine any patterns of confusion, 
% create a matrix that shows for each social message label (y-axis) the number of "yes" responses to each of the different 
% facial expression types (x-axis). If subjects tend to be accurate, the diagonal squares will have high values. 
% Display your results as a color-coded matrix using a different color map from Q1. 
% To compare data across cultures, use a 1 X 2 (cultures) subplot. 
% Include all necessary features (e.g., title, axes titles, colorbar). 
% Save this figure as a .jpg file in the Results Figures folder (3 POINTS)


% reshape matrix into more useable format 
data_rejigged = reshape(permute(all_data, [1, 3, 2, 4]),[nr_trials .* nr_subjects, 4, nr_cultures]);

% prepare matrix for data storage
confusion_matrix = zeros(nr_social_messages, nr_social_messages, nr_cultures);


for culture = 1:nr_cultures
    
    % find all yes responses
    yes_responses = data_rejigged(data_rejigged(:, 3, culture) == 1, :, culture);
    nr_trials_all = size(yes_responses, 1);
    
    % organise each trial into matrix
    for trial = 1:nr_trials_all
        
        facial_expression_type = yes_responses(trial, 1); % stimulus category
        word = yes_responses(trial, 2); % word presented with stimulus        
        confusion_matrix(word, facial_expression_type, culture) = confusion_matrix(word, facial_expression_type, culture) + 1;
    end    
end 


%% display figure 

max_yes = max(confusion_matrix(:));

figure, 
for culture = 1:nr_cultures

    subplot(1, 2, culture)
    imagesc(confusion_matrix(:, :, culture), [0 max_yes])
    
    axis image
    set(gca, 'YTick', 1:nr_social_messages, 'YTickLabel', social_messages)
    set(gca, 'XTick', 1:nr_social_messages, 'XTickLabel', social_messages)
    ylabel('Words')
    xlabel('Facial expression')
    
    title(cultures{culture})
    colorbar
end

% change colormap, print and save figure 
colormap hot 
cd([od '/' new_folder_name])
print('Confusion_matrices', '-djpeg')


%% Q5. Understanding the source of low d-prime values (a). 
% D-prime values can reflect 1 of four response patterns

% (a) high hits, low false alarms (facial expressions transmit the intended social message to the subject), 
% (b) low hits, high false alarms (the facial expressions transmit another social message to the subject than the one intended; they are confused with another message), 
% (c) high hits, high false alarms (the facial expressions transmit more than one social message to the subject; they are ambiguous), 
% (d) low hits, low false alarms (the facial expressions do not transmit any social message to the subject). 
% For each facial expression type and culture separately, create a figure that shows how each subject's responses distribute across 
% these four different response categories
% Include all necessary features (e.g., title, axes titles, colorbar). (3 POINTS)

culture_color = {'r', 'b'};

figure,
for facial_expression_type = 1:nr_social_messages
    for culture = 1:nr_cultures
        
        % extract hits and false alarm rate
        hits_prop = hits_fa_proportion(:, 1, facial_expression_type, culture);
        fa_prop = hits_fa_proportion(:, 2, facial_expression_type, culture);
        
        % plot colored points per culture
        subplot(2, 2, facial_expression_type)
        plot(fa_prop, hits_prop, 'o', 'MarkerFaceColor', culture_color{culture}, 'MarkerEdgeColor', 'w')
        hold on 
        
        axis([0 1 0 1])
        xlabel('False Alarm Rate')
        ylabel('Hit Rate')
    end
    
    title(social_messages{facial_expression_type})
    
    % add a dashed line for the four quadrants
    hold on 
    x = 0:.1:1;
    y = repmat(.5, [1 length(x)]);
    plot(x, y, 'LineStyle', '--', 'Color', 'w')
    
    hold on 
    y = x;
    x = repmat(.5, [1 length(x)]);
    plot(x, y, 'LineStyle', '--', 'Color', 'w')
    
end
legend(cultures, 'Location', 'SouthEast')

%% Q6. Understanding the source of d-prime values (b). 
% A subject's d-prime value could also be diminished or increased by specific facial expressions that omit or include 
% specific face movements (Action Units - AUs) that are necessary for accurate recognition. To examine this, show how 
% d-prime varies as a function of specific AU presence or absence in the facial expressions. 

% To do this, compute for each facial expression type and culture separately the d-prime value of each individual 
% facial expression by pooling the responses of the subjects. Then for each individual AU, compute separate average 
% d-prime values for when the AU is present and when it is absent, and compute the difference between them. 
% Plot your results using the y-axis for the AUs and the x-axis for the
% difference in d-prime values. Using subplotting and color-coding as you see fit. 
% Add a white line demarcating the point of zero difference. 
% Include all necessary features (e.g., title, axes titles, colorbar). 
% Finally, for each facial expression type and culture separately, rank the AUs according to their diagnosticity 
% - that is the magnitude of the increase in d-prime when they are present
% - and show a list of their AU names in the command line 
%(4 POINTS)

cd(od)
load facial_expression_models
load long_names42
[~, nr_aus, ~, nr_models] = size(data_EA_20sj);

dp_facial_expressions = zeros(nr_models, nr_social_messages, nr_cultures);
cultures_short = {'EA', 'WC'};

figure, 
for culture = 1:nr_cultures
    for facial_expression_type = 1:nr_social_messages
        for model = 1:nr_models
            
            % extract data for facial expression type and specific model
            loc_data = data_rejigged(:, 1, culture) == facial_expression_type & data_rejigged(:, 4, culture) == model;
            data = data_rejigged(loc_data, 2:3, culture);
            
            % number of matches and mismatches
            nr_matches = sum(data(:, 1) == facial_expression_type);
            nr_mismatches = sum(data(:, 1) ~= facial_expression_type);
            
            % hits and false alarms 
            hits = sum(data(:, 1) == facial_expression_type & data(:, 2) == 1);
            false_alarms = sum(data(:, 1) ~= facial_expression_type & data(:, 2) == 1);
            
            % compute proportions
            hits_proportion = hits ./ nr_matches;
            FA_proportion = false_alarms ./ nr_mismatches;
            
            % compute d-prime per model with downloaded function
            dp_facial_expressions(model, facial_expression_type, culture) = dprime_simple(hits_proportion, FA_proportion);
        end
        
        
        % compute average d-prime across all AUs when present versus absent
        mn_dp_au = zeros(nr_aus, 2);
        
        for au = 1:nr_aus
            eval(['loc_model_au = logical(squeeze(data_' cultures_short{culture} '_20sj(1, au, facial_expression_type, :)));'])
            
            % absent
            mn_dp_au(au, 1) = median(dp_facial_expressions(~loc_model_au, facial_expression_type, culture));
            
            % present
            mn_dp_au(au, 2) = median(dp_facial_expressions(loc_model_au, facial_expression_type, culture)); 
        end
        
        % positive - present is higher than absent 
        % negative - absent is higher than present 
        difference_dp = diff(mn_dp_au, 1, 2);
        
        subplot(2, 2, facial_expression_type)
        y = difference_dp;
        x = 1:nr_aus;
        plot(y, x, 'o', 'MarkerFaceColor', culture_color{culture}, 'MarkerEdgeColor', 'w')
        hold on
        set(gca, 'YTick', 1:nr_aus, 'YTickLabel', long_names42)
        title(social_messages(facial_expression_type))
        
        if facial_expression_type == 2
            legend(cultures, 'Location', 'NorthWest')
        end
        
        % add dashed line
        hold on
        x = zeros(1, nr_aus);
        y = 1:nr_aus;
        plot(x, y, 'LineStyle', '--', 'Color', 'w')
        
        % rank AUs
        [sorted, ind] = sort(difference_dp, 'descend');
        long_names42(ind)
    end
end
 
%% Q7. Sharing your data. 
% You need to send all of your data to a collaborator in another lab. 
% Since they do not use MATLAB, they need all the data prepared and saved as .xls files. 
% Prepare and save all of your data as .txt files using a sensible arrangement. (3 POINTS)           

% 1.	facial_expression_models.mat (per culture)
% (a)	Temporal parameters: (1): AU is present ? 1, absent ? 0, (2) AU peak amplitude, (3) AU peak latency, (4) AU onset latency, (5), AU offset latency, (6) AU deceleration, (7) AU acceleration. Values for 2 -7 have been normalized to the range 0-1.
% (b)	AUs: Corresponding AU names can be found in long_names42.mat
% (c)	Social messages: 1 ? thinking, 2 ? interested, 3 ? bored, and 4 ? confused
% (d)	Individual facial expressions. Location corresponds with the numbers in model_num (see 3d below)



new_dir = [od '/Files_for_my_collaborator'];
mkdir(new_dir)
cd(new_dir)

temporal_parameters = {'on_off', 'peak_amplitude', 'peak_latency', 'onset_latency', 'offset_latency', 'deceleration', 'acceleration'};
response_labels = {'yes', 'no'};


for culture = 1:nr_cultures
    for facial_expression_type = 1:nr_social_messages
        
        % 1. stimulus information for each observer in each culture
        for model = 1:nr_models
            
            % extract the facial expression model
            eval(['AU_pattern = data_' cultures_short{culture} '_20sj(:, :, facial_expression_type, model)''' ';'])
            stimulus_information = array2table(AU_pattern, 'RowNames', long_names42, 'VariableNames', temporal_parameters);
            
            fname = [cultures{culture} '_stimulus_information_model_' num2str(model) '.txt'];
            writetable(stimulus_information, fname); % saves as .txt file
        end
        
        
        % 2. Behavioural responses
        % prepare cell array space 
        
        face_on_trial = cell(nr_trials, 1);
        word_on_trial = face_on_trial;
        response_on_trial = face_on_trial;
        
        
        for subjects = 1:nr_subjects
            
            % model numbers 
            model_numbers = all_data(:, 4, subject, culture);
            
            % convert nuumbers to words where necessary 
            stimulus_category = all_data(:, 1, subject, culture);
            word_presented = all_data(:, 2, subject, culture);
            response_type = all_data(:, 3, subject, culture);
            
            for stim_cat = 1:nr_social_messages
                
                % facial expression 
                loc = find(stimulus_category == stim_cat)';
                face_on_trial(loc) = repmat(social_messages(stim_cat), size(loc));
                
                % word on trial
                loc = find(word_presented == stim_cat)';
                word_on_trial(loc) = repmat(social_messages(stim_cat), size(loc));
            end
            
            % yes/no response 
            for resp_type = 1:length(response_labels)
                loc = find(response_type == resp_type);
                response_on_trial(loc) = repmat(response_labels(resp_type), size(loc));
            end 
          
            % write to file 
            behavioural_data = table(model_numbers, face_on_trial, word_on_trial, response_on_trial);
            fname = [cultures{culture} '_behavioural_data_subject_' num2str(subjects) '.txt'];
            writetable(behavioural_data, fname); % saves as .txt file
        end
    end
end

%% End of exam