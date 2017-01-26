# ExtendingRipleysK

This includes the functions and scripts used to extend Ripley's K- function to quantify aggregation in 2-D grayscale images. 

Corresponding publication:

Extending Ripley's K-Function to Quantify Aggregation in 2-D Grayscale Images.
Amgad M, Itoh A, Tsui MM. PLoS One. 2015 Dec 4;10(12):e0144404. doi: 10.1371/journal.pone.0144404. eCollection 2015.

If you are unfamiliar with computer programming, we advise you to read through the following instructions to help you calculate K~ and the critical quantiles (Q01 and Q99) for your images at the desired range of radii. Beside the instructions below, you will find a series of seven screenshots showing, in a step-by-step manner, how to use the main script to batch process images in the input directory.

1- Open MATLAB
2- In the command window, type the following:
	cd('CodeDirectoryHere')
Of course, instead of 'CodeDirectoryHere', you would input the location of the supplementary file (S2), which contains all the functions and scripts accompanying the publication.
3- Strike the return button ("Enter"), and you will notice that the Current Folder at the left side of the screen now contains all the functions included in the S2 file.
4- Open the script called "CalculateKProfile.m"
5- Click "Run" (at the top screen, under the "Editor" bar).
6- A pop-up window will appear, asking for input parameters. Make sure you specify the input directory (containing your images) and the output directory where the K~ plots and the Excel file containing your results will be saved. If you do NOT want zero pixels to be considered as part of the field of view, type "1" in the corresponding field. Otherwise, leave as zero. 
7- Click "OK", and observe the Comman Window to monitor the progress of the calculations.
8- Once finished, a double arrow will show at the bottom of the text in the Command Window.
9- Go to the output directory and check your results.

Kindly feel free to contact the corresponding authors if you need any help or you have any advices or comments.
