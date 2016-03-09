Readme 
This is an image registration Matlab program. It contains two main folders, the 2D_multiresolution, and the 3D_multiresolution. This code works with most of the images.  Other than the given examples, you are free to try it with any images by saving the images in the corresponding folder and assign them as I0 (moving image) and I1 (fixed image).  Find out more detail on how to use this code by reading the Description_2D.docx and Description_3D.docx. 
Part of the code used the parfor loop which accelerates the program by performing parallel computing. However, it would require the user to have the license for Parallel Computing Toolbox. Also, the code can further speed up by using Mex function. You can generate the relevant Mex functions for this program by running the scripts make_mex_files2D.m (in the folder 2D_multiresolution) and make_mex_files3D.m (In the folder 3D image registration). The Mex function generated using these scripts are ready to use; you can then run the main file in the corresponding folder directly and enjoy the accelerated program. Short description about important variables used in the code can be find in variableList2D.docx and variableList3D.docx 


  
 



