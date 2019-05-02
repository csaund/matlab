% You work for the Glasgow City Council and are responsible for analysing data that 
% has been recorded by a CCTV camera in the city centre. This CCTV camera has filmed 
% a short section of Sauchiehall Street. The council wants to estimate how many 
% people pass through this section to understand how busy this part of the city
% centre is. Using this knowledge, the council will then decide which parts of 
% the city centre should be invested in for pedestrian use.
% You have been given a CCTV recording of 1 minute (recorded at 25 frames per second)
% to analyse. This CCTV-recorded data has already been processed by a software that 
% detects people. The output of this software is a series of single frames (i.e., 
% static images) where positive integers in each frame indicate that a human has 
% been detected. Within a given frame, each detected human has been assigned a 
% unique positive integer (e.g., 4). However, the human detection software cannot 
% connect consecutive frames. This means that the same person can be assigned a 
% different number across different frames (e.g., 1 in one frame, but 5 in the next 
% frame), which makes counting the number of people who have passed through this 
% section problematic.
% Your task is to resolve this problem by applying a simple heuristic (details 
% provided below) so that each person detected keeps the same identifier across
% frames.
% Data files
% You have been given a folder which contains the following:
% 1. imagesList.mat This file contains a list of the image filenames recorded by
% the CCTV.
% 2. imageFiles,afoldercontaining.pngfiles.Thesearetheimagefilesrecorded
% by the CCTV. The second number in the file names indicates the temporal
% order of the images starting at C0039_100000_INDS.png, i.e frame number 100000
% is the first image).

%% 
%% setup
clear all

cd('/Users/carolynsaund/Desktop/exam/matlab_exam_2019/Data')
od = pwd;



%% 
% Q1 (2 points)
% First, make sure that you have received all of the relevant image files. 
% Check whether all image files in the list are included in the folder of image files. 
% How many images are there? Are any frames missing? Display the results in the command
% window.
% 1 minute = 60s x 25 fps = 1500 images 

%% 

% Q2 (3 points)
% Next, get an impression of what the images look like. To do this, randomly select 
% two frames with consecutive numbers AND with people detected in them. Use any image 
% display function to present these images in a single figure and positioned side-by- 
% side in a 1 x 2 subplot.

%% 

% Q3 (3 points)
% Now, get an overview of how many people are detected across the images. 
% To do this, count the number of different people detected in each frame and
% display the results as a histogram.

%% 

% Q4 (4 points)
% Your task is to estimate how many people walked through this section of the city 
% centre during the recording. However, in each frame, each detected person has been 
% assigned a unique number, which makes this estimation problematic. Your job now is
% to infer whether, in consecutive frames, a given identifier represents the same 
% person or a different person.
% To address this problem, calculate the amount of overlap between pixels that 
% identify a given person in one frame and the pixels that identify a given person
% in the previous frame (i.e., how many pixels of person with ID number 4 in frame 2
% overlap with the person with ID number 5 in frame 1?). If there is high overlap 
% between detected persons, this likely indicates that these are in fact the same 
% person. Divide this amount by the number of pixels that the respective person 
% occupies in the current frame (i.e., calculate the percentage of pixels that
% overlap). Calculate this
%  Make sure the proportions of the plots retain the original
%  proportions of the images. Use the ?hot? colormap. Add the time point of the frames
% as the title of each image (e.g., time point 1, time point 3).
% value for all pairs of consecutive frames ([1, 2], [2, 3], [3, 4] ...) and for all
% detected persons. Display the pooled results as a histogram using 100 bins. Crop 
% the figures on the y-axis so that you can better see the shape of the distribution
% and to ignore most of the entries for the lowest bin.

%% 

% Q5 (6 points)
% Next, you must decide whether each person detected in each frame corresponds to the
% same person in the neighbouring frame, or if it?s a different person. To do this, 
% choose a cut-off value that roughly divides the histogram created in Q5 in two 
% (don?t worry ? your choice of cut-off value will not be marked but do make sure 
% that the cut- off is sensible). Now, for each pair of consecutive frames and for 
% each detected person in each frame, identify any persons in the preceding frame 
% whose level of pixel overlap exceeds that threshold. For any such persons, assign 
% the identity of the corresponding person in the previous frame to the person in 
% the current frame.
% Next, assume that a person cannot walk through this section of the city centre in 
% less than 3 seconds. Using this rule, treat all identified persons that appear to
% pass through this section in less than 3 seconds as noise and remove them.
% Display your results in a figure as a matrix where the x-axis corresponds to time 
% in seconds and the y-axis corresponds to unique persons. The appearance of each 
% person across time (i.e., across frames) will be represented as a 1 in the matrix,
% all other values as 0. Use the colormap ?gray?. Add to the title the final inferred 
% number of persons that have passed through this section of the city centre, the 
% mean time needed to pass through this section and its standard deviation.
% Finally, display a line plot showing the total number of persons that passed 
% through this section of the city centre at each time point. Use a thick red line.
% Use solid green circles to mark each time point. Add meaningful x- and y-axis 
% labels and a title.

%% 

% Q6 (2 points)
% To share your results with your colleague, create a folder that contains three 
% .mat- files: 
% (1) the cut-off value you chose in Q5; 
% (2) a structure containing 
% (a) each raw frame, 
% (b) the corresponding matrix containing the numerical identifiers of 
% each person detected in the frame after you have applied the heuristic
% in Q5 (i.e., the processed image), 
% (c) a value indicating the number of people detected in that frame, and 
% (d) a character array (i.e., string) indicating the time point of that 
% frame (e.g., ?Time point 10?); and 
% (3) a 3D matrix containing the processed images where the third dimension is time.

%% 

% Q7 (2 points)
% Save all created figures as coloured pdf files in another folder. Set the sizes of the figures so that the figure contents are clearly visible and that all text items (axis labels, titles, etc) are easy to read.
